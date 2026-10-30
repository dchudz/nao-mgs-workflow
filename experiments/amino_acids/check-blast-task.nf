include { BLASTN as ORIGINAL } from './baseline-blast'
include { BLASTN as DEFAULT } from '../../modules/local/blast'
include { BLASTN as SENSITIVE } from '../../modules/local/blast'
workflow {
    input = channel.of(tuple('fixture',file(params.query,checkIfExists:true)))
    settings = [blast_perc_id:60,blast_qcov_hsp_perc:30,db_download_timeout:1200]
    ORIGINAL(input,params.database,settings)
    DEFAULT(input,params.database,settings)
    SENSITIVE(input,params.database,settings+[blast_task:'blastn'])
}
