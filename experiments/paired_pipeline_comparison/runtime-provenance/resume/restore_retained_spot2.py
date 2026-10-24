
import subprocess,json,pathlib,datetime
rootdev='/dev/disk/by-id/nvme-Amazon_Elastic_Block_Store_vol05ca1445d265f5325-part1'
scratchdev='/dev/disk/by-id/nvme-Amazon_Elastic_Block_Store_vol03f6a458f0812c7b4'
old=pathlib.Path('/mnt/retained-spot2')
def run(*args):subprocess.run(list(args),check=True)
assert not subprocess.check_output(['sudo','docker','ps','-q'],text=True).strip()
run('sudo','systemctl','stop','docker')
run('sudo','mkdir','-p',str(old));run('sudo','mount','-o','nouuid',rootdev,str(old))
assert (old/'home/ec2-user/paired-study-20261024/source-commit.txt').read_text().strip()=='c77cdde08b03834a9c1f54f3d80b0c9668bf9a83'
run('sudo','mount','--bind',str(old/'home/ec2-user'),'/home/ec2-user')
run('sudo','mount','--bind',str(old/'var/lib/docker'),'/var/lib/docker')
run('sudo','mkdir','-p','/mnt/paired-study-inputs');run('sudo','mount',scratchdev,'/mnt/paired-study-inputs')
# Persist mounts by unique EBS device identity; do not use the duplicated root UUID.
entries=[rootdev+' /mnt/retained-spot2 xfs defaults,nouuid 0 0',str(old/'home/ec2-user')+' /home/ec2-user none bind,x-systemd.requires-mounts-for=/mnt/retained-spot2 0 0',str(old/'var/lib/docker')+' /var/lib/docker none bind,x-systemd.requires-mounts-for=/mnt/retained-spot2 0 0',scratchdev+' /mnt/paired-study-inputs xfs defaults 0 0']
existing=pathlib.Path('/etc/fstab').read_text();assert '/mnt/retained-spot2' not in existing
p=subprocess.run(['sudo','tee','-a','/etc/fstab'],input='\n'+'\n'.join(entries)+'\n',text=True,capture_output=True,check=True)
run('sudo','systemctl','daemon-reload');run('sudo','systemctl','start','docker')
images={image:subprocess.check_output(['sudo','docker','image','inspect',image,'--format','{{.Id}}'],text=True).strip() for image in ['pipeline-e2e-rust:local','public.ecr.aws/q0n1c7g8/nao-mgs-workflow/pandas:91e85c3fd8e9026c']}
assert images['pipeline-e2e-rust:local']=='sha256:d1b77ec75262ad7dc27dce983b4294782c84d986529aabe335c94d8f9aab7b12'
assert not subprocess.check_output(['sudo','docker','ps','-q'],text=True).strip()
print(json.dumps({'restored_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'root_data_volume':'vol-05ca1445d265f5325','scratch_volume':'vol-03f6a458f0812c7b4','images':images,'mounts':subprocess.check_output(['findmnt','-J','/home/ec2-user','/var/lib/docker','/mnt/paired-study-inputs'],text=True)}))
