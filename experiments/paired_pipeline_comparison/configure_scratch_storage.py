"""Mount a verified new, empty AWS volume for temporary study inputs."""
import json
import os
import subprocess
import sys
from pathlib import Path

worker, volume_id = sys.argv[1:]
serial = volume_id.replace('-', '')
mount = Path('/mnt/paired-study-inputs')
items = json.loads(subprocess.check_output(['lsblk','-J','-d','-o','NAME,SERIAL,FSTYPE,SIZE'],text=True))['blockdevices']
matching = [x for x in items if x['serial'] == serial]
assert len(matching) == 1
x = matching[0]
assert x['name'] != 'nvme0n1' and x['size'] == '512G' and not x['fstype']
dev = '/dev/'+x['name']
probe = subprocess.run(['sudo','blkid','-p',dev],capture_output=True,text=True)
assert probe.returncode == 2 and not probe.stdout
subprocess.run(['sudo','mkfs.xfs',dev],check=True,stdout=subprocess.DEVNULL)
subprocess.run(['sudo','mkdir','-p',str(mount)],check=True)
assert not list(mount.iterdir())
subprocess.run(['sudo','mount',dev,str(mount)],check=True)
subprocess.run(['sudo','chown',f'{os.getuid()}:{os.getgid()}',str(mount)],check=True)
uuid = subprocess.check_output(['sudo','blkid','-s','UUID','-o','value',dev],text=True).strip()
fstab = subprocess.check_output(['sudo','cat','/etc/fstab'],text=True)
assert str(mount) not in fstab
line = f'UUID={uuid} {mount} xfs defaults,nofail,x-systemd.device-timeout=30 0 2\n'
subprocess.run(['sudo','tee','-a','/etc/fstab'],input=line,text=True,check=True,stdout=subprocess.DEVNULL)
result = dict(worker=worker,volume_id=volume_id,device=dev,filesystem_uuid=uuid,mount=str(mount),
              empty_volume_verified_before_format=True,
              purpose='Temporary source acquisition and sampling only; reference and pipeline output disks unchanged.')
Path('/home/ec2-user/paired-study-20261024/scratch-storage.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
