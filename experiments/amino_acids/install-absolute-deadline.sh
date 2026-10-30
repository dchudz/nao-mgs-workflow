#!/bin/bash
set -euo pipefail
sudo tee /etc/systemd/system/aa-hard-deadline.service >/dev/null <<'UNIT'
[Unit]
Description=Terminate authorized experiment at absolute deadline
[Service]
Type=oneshot
ExecStart=/sbin/shutdown -h now
UNIT
sudo tee /etc/systemd/system/aa-hard-deadline.timer >/dev/null <<'UNIT'
[Unit]
Description=Absolute experiment cleanup deadline
[Timer]
OnCalendar=2026-10-22 15:12:00 UTC
Persistent=true
Unit=aa-hard-deadline.service
[Install]
WantedBy=timers.target
UNIT
sudo systemctl daemon-reload
sudo systemctl enable --now aa-hard-deadline.timer
systemctl list-timers aa-hard-deadline.timer --no-pager
