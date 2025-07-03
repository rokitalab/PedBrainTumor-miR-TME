class: CommandLineTool
cwlVersion: v1.2
id: extract_3utr_genocde
doc: |-
  Extract 3' UTR sequence from Gencode.

requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.1
- class: ResourceRequirement
  ramMin: $(inputs.ram * 1000)
  coresMin: $(inputs.cores)
- class: InitialWorkDirRequirement
  listing:
  - entryname: 00-extract-3utr-gencode-v39.R
    writable: false
    entry:
      $include: ../00-extract-3utr-gencode-v39.R
baseCommand: [Rscript]
arguments:
- position: 1
  shellQuote: false
  valueFrom: |
    --output gencode.v39.3utr.fa
inputs:
  gencode_gtf: { type: File, inputBinding: { prefix: "--gencode", position: 1 }, doc: "Input gencode file" }
  ram: { type: 'int?', default: 8, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  array_dirs:
    type: 'File'
    outputBinding:
      glob: "gencode.v39.3utr.fa"
