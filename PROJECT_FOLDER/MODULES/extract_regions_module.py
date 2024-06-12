from gwf import AnonymousTarget

def get_observable_sequence(phased_file_dir, phased_file_prefix, outfile_dir):
    """
    FIXME
    """
    
    inputs = [f'{phased_file_dir}/{phased_file_prefix}_1.bed', f'{phased_file_dir}/{phased_file_prefix}_2.bed']
    outputs = [f'{outfile_dir}/temp_dir/split.DONE']
    options = {"walltime":"1:00:00","account":"sexChromosomes", "memory":"20gb", "cores": 16}
    spec = f'''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pythonHMM

    echo "Job ID: $SLURM_JOB_ID"
    
    mkdir -p {outfile_dir}/temp_dir 
    
    python SCRIPTS/extract_observable_sequence.py \
        --infiles {f'{phased_file_dir}/{phased_file_prefix}_1.bed'} {f'{phased_file_dir}/{phased_file_prefix}_2.bed'} \
        --outdir {outfile_dir}/temp_dir \
        --split_chromosomes
    
    touch {outfile_dir}/temp_dir/split.DONE
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def get_imprinted_regions(infile, outfile_dir, chrom, parameter_file):
    """
    FIXME
    """
    
    inputs = [f'{outfile_dir}/temp_dir/split.DONE']
    outputs = [f'{outfile_dir}/temp_dir/{chrom}_imprinting_HMM.DONE']
    options = {"walltime":"48:00:00","account":"sexChromosomes", "memory":"5gb", "cores": 4}
    spec = f'''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pythonHMM

    echo "Job ID: $SLURM_JOB_ID"
    
    python SCRIPTS/run_imprinting_HMM.py \
        --infile {infile} \
        --outdir {outfile_dir} \
        --parameters {parameter_file} \
        --max_iterations 20 \
        --min_iterations 5 \
        --stop_condition 1
    
    touch {outfile_dir}/temp_dir/{chrom}_imprinting_HMM.DONE
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def merge_regions(infile, outfile, regions, meta_file, extra_flags):
    """
    FIXME
    """
    
    inputs = [regions]
    outputs = [outfile]
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"264gb", "cores": 4}
    spec = f'''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pythonHMM

    echo "Job ID: $SLURM_JOB_ID"
    
    echo "python SCRIPTS/aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {regions} \
        --sample_meta {meta_file} \
        {extra_flags}"
    
    python SCRIPTS/aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {regions} \
        --sample_meta {meta_file} \
        {extra_flags}
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
