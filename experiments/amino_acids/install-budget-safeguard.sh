#!/bin/bash
set -euo pipefail
sudo tee /etc/systemd/system/aa-budget-safeguard.service >/dev/null <<'UNIT'
[Unit]
Description=Terminate experiment compute at conservative budget boundary
[Service]
Type=oneshot
ExecStart=/sbin/shutdown -h now
UNIT
sudo tee /etc/systemd/system/aa-budget-safeguard.timer >/dev/null <<'UNIT'
[Unit]
Description=Experiment budget safeguard
[Timer]
OnCalendar=2026-10-23 15:00:00 UTC
Persistent=true
Unit=aa-budget-safeguard.service
[Install]
WantedBy=timers.target
UNIT
sudo systemctl daemon-reload
sudo systemctl enable --now aa-budget-safeguard.timer
sudo systemctl stop aa-experiment-deadline.timer 2>/dev/null || true
sudo systemctl disable --now aa-hard-deadline.timer 2>/dev/null || true
systemctl list-timers --all --no-pager | grep -E 'aa-|experiment|shutdown'
