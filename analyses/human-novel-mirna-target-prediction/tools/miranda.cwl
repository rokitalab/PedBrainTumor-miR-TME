class: CommandLineTool
cwlVersion: v1.2
id: miranda
doc: |-
  Run miranda

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.1
- class: ResourceRequirement
  ramMin: $(inputs.ram * 1000)
  coresMin: $(inputs.cores)
baseCommand: [miranda]
arguments:
- position: 99
  prefix: ''
  shellQuote: false
  valueFrom: |
    results 1>&2
inputs:
  target_fasta : {type: File, inputBinding: {position: 1}, doc: "Fasta file with target miRNA queries"}
  mirna_fasta : {type: File, inputBinding: {position: 2}, doc: "Fasta file with sequence data"}
  ram: { type: 'int?', default: 32, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 8, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  mirand_output_file:
    type: 'File'
    outputBinding:
      glob: 'results/miranda_output.txt'
