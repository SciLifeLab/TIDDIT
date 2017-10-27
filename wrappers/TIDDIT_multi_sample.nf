//absolute path to the TIDDIT binary
TIDDIT_path="/home/jesper/TIDDIT/bin/TIDDIT"
//absolute path to the TIDDIT_combine.py script
combine_script_path="/home/jesper/TIDDIT/wrappers/TIDDIT_multi_sample.py"

//default output directory
params.working_dir="TIDDIT_multi_sample"

//variants are called at these limits
call_read_pairs=3
call_split_reads=2

//variants are kept if the pairs/split reads of a sample exceed these limits
report_read_pairs=6
report_split_reads=4

//overlap and distance to concider two variants the same
bnd_distance = 1000
overlap = 0.5


if(!file(TIDDIT_path).exists()) exit 1, "could not find the TIDDIT binary, make sure to set the TIDDIT_path variable"
if(!file(combine_script_path).exists()) exit 1, "could not dine the TIDDIT_multi_sample.py script, please set the combine_script_path variable"
params.bam=""
if(!params.bam) exit 1, "error no bam files, enter the bam files separated by , example --bam file1.bam,file2.bam"
params.output=""
if(!params.output) exit 1, "error no output file name, set the output parameter files separated by , --out output.vcf"

Channel.from( params.bam.splitCsv()).subscribe{if(!file(it).exists()) exit 1, "Missing bam:${it}"}

bam_files=Channel.from(params.bam.splitCsv()).map{
    line -> file(line)
} 

process TIDDIT {
    publishDir "${params.working_dir}", mode: 'copy', overwrite: true     
    tag { bam_file }
    
    cpus 1
        
    input:
    file(bam_file) from bam_files
    
    output:
    file "${bam_file.baseName}.vcf" into TIDDIT_vcf
    
    script:
    """
    ${TIDDIT_path} --sv -b ${bam_file} -p ${call_read_pairs} -r ${call_split_reads} -q 10 -o ${bam_file.baseName}
    rm *.tab
     """
}

process combine {
    publishDir "${params.working_dir}", mode: 'copy', overwrite: true
    
    cpus 1
        
    input:
    file vcf from TIDDIT_vcf.toList()
    
    output:
    file "${params.output}" into multi_sample_vcf
    
    script:
    """
    svdb --merge --bnd_distance ${bnd_distance} --overlap ${overlap} --vcf *.vcf > tmp.vcf
    python ${combine_script_path} tmp.vcf ${report_split_reads} ${report_read_pairs} > ${params.output}
    """
}
