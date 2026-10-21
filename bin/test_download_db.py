#!/usr/bin/env python3
"""
Unit tests for download-db.py

Run with: pytest bin/test_download_db.py
"""

import fcntl
import os
import stat
import subprocess
from collections.abc import Generator
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from download_db import (
    download_database,
    file_lock,
    get_cache_name,
    main,
    parse_source_path,
    sync_from_local,
    sync_from_s3,
)


class TestParseSourcePath:
    """Test source path parsing and normalization."""

    @pytest.mark.parametrize(
        "input_path,expected_path,expected_is_s3",
        [
            # S3 paths - various malformed inputs
            ("s3://bucket/path/db", "s3://bucket/path/db", True),
            ("s3:///bucket/path/db", "s3://bucket/path/db", True),
            ("s3:/bucket/path/db", "s3://bucket/path/db", True),
            ("s3:bucket/path/db", "s3://bucket/path/db", True),
            ("s3://bucket//path///db", "s3://bucket/path/db", True),
            # Local paths
            ("/path/to/db", "/path/to/db", False),
            ("/path//to///db", "/path/to/db", False),
            ("relative/path/db", "relative/path/db", False),
        ],
    )
    def test_parse_source_path(
        self, input_path: str, expected_path: str, expected_is_s3: bool
    ) -> None:
        """Test path normalization with various inputs."""
        normalized, is_s3 = parse_source_path(input_path)
        assert normalized == expected_path
        assert is_s3 == expected_is_s3


class TestGetCacheName:
    """Test cache name generation."""

    @pytest.mark.parametrize(
        "source_path,expected_name",
        [
            ("s3://bucket/path/kraken_db", "kraken_db_35f571f1"),
            ("s3://bucket/path/db", "db_56898535"),
            ("/local/path/bowtie2_idx", "bowtie2_idx_46a0d77f"),
        ],
    )
    def test_cache_name(self, source_path: str, expected_name: str) -> None:
        """Test that cache name is basename + first 8 chars of sha256 hash."""
        assert get_cache_name(source_path) == expected_name


class TestFileLock:
    """Test file locking context manager."""

    @pytest.mark.parametrize("lock_path", ["test.lock", "subdir/nested/test.lock"])
    def test_lock_creates_file_and_parents(
        self, tmp_path: Path, lock_path: str
    ) -> None:
        """Test that lock creates the lock file and any parent directories."""
        lock_file = tmp_path / lock_path
        with file_lock(lock_file):
            assert lock_file.exists()

    def test_lock_timeout(self, tmp_path: Path) -> None:
        """Test that timeout raises TimeoutError when lock is held by another."""
        lock_file = tmp_path / "test.lock"
        with open(lock_file, "w") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            with (
                pytest.raises(TimeoutError, match="Timed out waiting for lock"),
                file_lock(lock_file, timeout_seconds=0.2),
            ):
                pass
            # Lock is released when file is closed


class TestSyncFromS3:
    """Test isolated S3 configuration and database sync."""

    @pytest.mark.parametrize("custom_config", [False, True])
    @patch("download_db.subprocess.run")
    def test_sync_preserves_user_configuration(
        self,
        mock_run: MagicMock,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        custom_config: bool,
    ) -> None:
        """Transfer settings are isolated; existing config and environment survive."""
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.delenv("AWS_CONFIG_FILE", raising=False)
        original = tmp_path / ".aws" / "config"
        if custom_config:
            original = tmp_path / "custom.config"
            monkeypatch.setenv("AWS_CONFIG_FILE", "$HOME/custom.config")
        monkeypatch.setenv("AWS_PROFILE", "research")
        original.parent.mkdir(parents=True, exist_ok=True)
        contents = (
            b"[profile research]\nsource_profile = base\nregion = us-east-1\n"
            b"[profile base]\ncredential_process = example-credential-helper\n"
        )
        original.write_bytes(contents)
        original_env = os.environ.copy()
        configs = []

        def run(command: list[str], *, env: dict[str, str], **kwargs: object) -> None:
            config = Path(env["AWS_CONFIG_FILE"])
            configs.append(config)
            assert config != original
            assert config.read_bytes() == contents
            assert stat.S_IMODE(config.stat().st_mode) == 0o600
            assert env == {**original_env, "AWS_CONFIG_FILE": str(config)}

        mock_run.side_effect = run
        sync_from_s3("s3://bucket/db", tmp_path / "db")
        assert len(set(configs)) == 1
        assert not configs[0].exists()
        assert original.read_bytes() == contents
        assert os.environ == original_env
        settings = dict(
            call.args[0][3:]
            for call in mock_run.call_args_list
            if call.args[0][1] == "configure"
        )
        assert settings == {
            "s3.preferred_transfer_client": "crt",
            "s3.multipart_chunksize": "16MB",
        }

    @patch("download_db.subprocess.run")
    def test_configurations_do_not_share_files(
        self, mock_run: MagicMock, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Overlapping transfers have independent temporary configurations."""
        monkeypatch.setenv("AWS_CONFIG_FILE", str(tmp_path / "missing"))
        configs = []

        def run(command: list[str], *, env: dict[str, str], **kwargs: object) -> None:
            if command[:3] == ["aws", "s3", "sync"]:
                config = Path(env["AWS_CONFIG_FILE"])
                configs.append(config)
                if len(configs) == 1:
                    sync_from_s3("s3://bucket/second", tmp_path / "second")
                assert config.is_file()

        mock_run.side_effect = run
        sync_from_s3("s3://bucket/first", tmp_path / "first")
        assert len(configs) == 2
        assert configs[0] != configs[1]
        assert all(not config.exists() for config in configs)

    @patch("download_db.subprocess.run")
    def test_sync_from_s3_success(
        self, mock_run: MagicMock, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Test successful S3 sync."""
        monkeypatch.setenv("AWS_CONFIG_FILE", str(tmp_path / "missing"))
        sync_from_s3("s3://bucket/db", Path("/scratch/db"))
        env = mock_run.call_args_list[0].kwargs["env"]
        assert mock_run.call_args.args[0] == [
            "aws",
            "s3",
            "sync",
            "s3://bucket/db",
            "/scratch/db",
            "--delete",
        ]
        assert mock_run.call_args.kwargs["env"] == env
        assert mock_run.call_args.kwargs["check"]
        assert not Path(env["AWS_CONFIG_FILE"]).exists()

    @pytest.mark.parametrize("failing_command", ["configure", "s3"])
    @patch("download_db.subprocess.run")
    def test_sync_from_s3_failure(
        self,
        mock_run: MagicMock,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        failing_command: str,
    ) -> None:
        """Configuration or download failures propagate and clean up the config."""
        monkeypatch.setenv("AWS_CONFIG_FILE", str(tmp_path / "missing"))

        def run(command: list[str], **kwargs: object) -> None:
            if command[1] == failing_command:
                raise subprocess.CalledProcessError(1, command)

        mock_run.side_effect = run
        with pytest.raises(subprocess.CalledProcessError):
            sync_from_s3("s3://bucket/db", Path("/scratch/db"))
        assert not Path(mock_run.call_args.kwargs["env"]["AWS_CONFIG_FILE"]).exists()


class TestSyncFromLocal:
    """Test local rsync function."""

    @patch("download_db.subprocess.run")
    def test_sync_from_local_success(self, mock_run: MagicMock) -> None:
        """Test successful local sync."""
        sync_from_local("/source/db", Path("/scratch/db"))
        mock_run.assert_called_once_with(
            ["rsync", "-a", "--delete", "/source/db/", "/scratch/db/"],
            check=True,
            capture_output=True,
            text=True,
        )

    @patch("download_db.subprocess.run")
    def test_sync_from_local_failure(self, mock_run: MagicMock) -> None:
        """Test local sync failure raises CalledProcessError."""
        mock_run.side_effect = subprocess.CalledProcessError(1, "rsync")
        with pytest.raises(subprocess.CalledProcessError):
            sync_from_local("/source/db", Path("/scratch/db"))


class TestDownloadDatabase:
    """Test main download_database function."""

    @pytest.fixture
    def mock_file_lock(self) -> Generator[MagicMock, None, None]:
        with patch("download_db.file_lock") as mock:
            mock.return_value.__enter__ = MagicMock()
            mock.return_value.__exit__ = MagicMock(return_value=False)
            yield mock

    @pytest.mark.parametrize(
        "source_path,expected_basename,sync_func",
        [
            ("s3://bucket/kraken_db", "kraken_db_", "download_db.sync_from_s3"),
            ("/source/bowtie2_db", "bowtie2_db_", "download_db.sync_from_local"),
            (
                "s3:///bucket//path///db",
                "db_",
                "download_db.sync_from_s3",
            ),  # tests normalization
        ],
    )
    def test_download_database(
        self,
        mock_file_lock: MagicMock,
        tmp_path: Path,
        source_path: str,
        expected_basename: str,
        sync_func: str,
    ) -> None:
        """Test downloading from S3 and local paths."""
        with patch(sync_func) as mock_sync:
            result = download_database(source_path, scratch_dir=tmp_path)
            assert result.parent == tmp_path
            assert result.name.startswith(expected_basename)
            mock_sync.assert_called_once()

    def test_download_database_creates_scratch_dir(
        self, mock_file_lock: MagicMock, tmp_path: Path
    ) -> None:
        """Test that scratch directory is created if it doesn't exist."""
        _ = mock_file_lock  # fixture needed for side effect
        with patch("download_db.sync_from_s3"):
            scratch = tmp_path / "new_scratch"
            assert not scratch.exists()
            download_database("s3://bucket/db", scratch_dir=scratch)
            assert scratch.exists()

    def test_download_database_passes_timeout(self, tmp_path: Path) -> None:
        """Test that timeout is passed to file_lock."""
        with (
            patch("download_db.file_lock") as mock_lock,
            patch("download_db.sync_from_s3"),
        ):
            mock_lock.return_value.__enter__ = MagicMock()
            mock_lock.return_value.__exit__ = MagicMock(return_value=False)
            download_database(
                "s3://bucket/db", timeout_seconds=600, scratch_dir=tmp_path
            )
            assert mock_lock.call_args[0][1] == 600


class TestMain:
    """Test main entry point."""

    @pytest.mark.parametrize(
        "argv,expected_args",
        [
            (["download_db.py", "s3://bucket/db", "300"], ("s3://bucket/db", 300)),
            (["download_db.py", "/local/path/db"], ("/local/path/db", None)),
        ],
    )
    @patch("download_db.download_database")
    def test_main(
        self,
        mock_download: MagicMock,
        argv: list[str],
        expected_args: tuple[str, int | None],
    ) -> None:
        """Test that main parses args and calls download_database."""
        mock_download.return_value = Path("/scratch/db_abc123")
        with patch("sys.argv", argv):
            main()
            mock_download.assert_called_once_with(*expected_args)

    @patch("download_db.download_database")
    def test_main_prints_path(
        self, mock_download: MagicMock, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Test that main prints the local path to stdout."""
        mock_download.return_value = Path("/scratch/kraken_db_abc12345")
        with patch("sys.argv", ["download_db.py", "s3://bucket/db"]):
            main()
        captured = capsys.readouterr()
        assert captured.out.strip() == "/scratch/kraken_db_abc12345"
