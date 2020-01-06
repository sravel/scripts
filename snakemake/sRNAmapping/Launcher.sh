#!/bin/bash
#rm snakejob.*
module purge
module load system/python/3.7.3

export DRMAA_LIBRARY_PATH=/SGE/8.1.8/lib/lx-amd64/libdrmaa.so
cluster_config="/work/chabannes/TEST/cluster_config.yaml"
datas_config="/work/chabannes/TEST/config.yaml"
scratch_dir="/tmp"

# produit le graph du pipeline
snakemake -s sRNA_pipeline.snake --configfile ${datas_config} --rulegraph  | dot -Tpdf > schema_pipeline_global.pdf



snakemake --shadow-prefix ${scratch_dir} --latency-wait 5184000 -s sRNA_pipeline.snake --jobs 100 --drmaa "{cluster.queue} {cluster.export_env} {cluster.cwd} {cluster.mem} {cluster.n_cpu}{threads} {cluster.logerror}{params.errorLog} {cluster.log}{params.outputLog} " --cluster-config ${cluster_config} --jobname "{rulename}.{name}.{jobid}" --configfile ${datas_config}


#snakemake --latency-wait 5184000 -s sRNA_pipeline.snake --jobs 100 --cluster "qsub -e {params.errorLog}.${JOBID} -o {params.outputLog}.${JOBID} -q {params.queue} -cwd -V -pe parallel_smp {threads} -l mem_free={params.l_mem_free}"

snakemake -s sRNA_pipeline.snake --configfile ${datas_config} --filegraph | dot -Tpdf > schema_pipeline_files.pdf
snakemake -s sRNA_pipeline.snake --configfile ${datas_config} --dag | dot -Tpdf > schema_pipeline_samples.pdf
snakemake -s sRNA_pipeline.snake --report REPORT.html

# add to clean files
#snakemake -s sRNA_pipeline.snake clean
