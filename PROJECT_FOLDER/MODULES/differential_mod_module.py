
from gwf import AnonymousTarget

def run_differential(sample_path_list, sample_flags, cpg_bedfile, log_file, prefix, outfile_dir, reference):
    """
    FIXME
    """
    inputs = sample_path_list
    outputs = [f'{outfile_dir}/{prefix}.bed'] ### Nok ikke rigtigt
    options = {"walltime":"24:00:00","account":"sexChromosomes", "memory":"200gb", "cores": 32}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit dmr multi {samples} --log-filepath {log} --prefix {prefix} --force --threads 32 --out-dir {out_dir} --ref {ref} --regions-bed {cpgs} --base C"
    
    modkit dmr multi \
        {samples} \
        --log-filepath {log} \
        --prefix {prefix} \
        --force \
        --threads 32 \
        --out-dir {out_dir} \
        --ref {ref} \
        --regions-bed {cpgs} \
        --base C
    
    '''.format(bed = " ".join(sample_path_list), samples = sample_flags, log = log_file, prefix = prefix, out_dir = outfile_dir, ref = reference, cpgs = cpg_bedfile)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


#nanopolish call-methylation [OPTIONS] --reads reads.fastq --bam alignments.bam --genome genome.fa --methylation cpg -v --threads 32 --batchsize 512
### FIXME: create metH5 files... return tsv pls

# $ meth5 create_m5 --input_paths INPUT_PATH1 [INPUT_PATH2 ...] --output_file OUTPUT_FILE.m5
# $ meth annotate_reads --m5file M5FILE.m5 --read_groups_key READ_GROUPS_KEY --read_group_file READ_GROUP_FILE

# $ meth5 convert -i INPUT_FILE1.bam INPUT_FILE2.bam INPUT_FILEN.bam -o OUTPUT_FILE.m5

def run_cpg_finder(bed_outfile, tsv_outfile, reference):
    """
    FIXME
    I'm using the default settings for:  [-m MERGE_GAP] [-w MIN_WIN_LEN] [-c MIN_CG_FREQ] [-r MIN_OBS_CG_RATIO]
    It can be in bed or tsv format: [-b OUTPUT_BED_FN] [-t OUTPUT_TSV_FN]
    """
    inputs = [reference]
    outputs = ["CpG_finder_test.DONE"] ### Nok ikke rigtigt
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"100gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONTmod

    echo "Job ID: $SLURM_JOB_ID"
    echo ""

    
    pycometh CGI_Finder -f {ref} -b {bed} -t {tsv} -v -p
    
    touch CpG_finder_test.DONE
    
    '''.format(bed = bed_outfile, tsv = tsv_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_meth_segmentation(bed_outfile, tsv_outfile, reference):
    """
    FIXME
    I'm using the default settings for:  [-m MERGE_GAP] [-w MIN_WIN_LEN] [-c MIN_CG_FREQ] [-r MIN_OBS_CG_RATIO]
    It can be in bed or tsv format: [-b OUTPUT_BED_FN] [-t OUTPUT_TSV_FN]
    """
    inputs = sample_path_list
    outputs = [f'{outfile_dir}/{prefix}.bed'] ### Nok ikke rigtigt
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"100gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONTmod

    echo "Job ID: $SLURM_JOB_ID"
    echo ""

    
    pycometh Meth_Seg -i {group_MetH5s} -r haplotype -c {chrom} -t {tsv} -r {group_key} -p 32 --print_diff_met
    
    touch Meth_seg_test.DONE
    
    '''.format(bed = bed_outfile, tsv = tsv_outfile, group_key = read_group_keys, MetH5 = MetH5_files, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



