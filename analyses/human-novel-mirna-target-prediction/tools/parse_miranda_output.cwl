class: CommandLineTool
cwlVersion: v1.2
id: parse_miranda_outputs
doc: |-
  Parse Miranda output file.

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
  - entryname: 00-unzip-and-sort.R
    writable: false
    entry:
      $include: ../03-parse-miranda-output.py
baseCommand: [python]
arguments:
- position: 99
  prefix: ''
  shellQuote: false
  valueFrom: |
    1>&2
inputs:
  input_miranda_results: { type: File, inputBinding: { prefix: "--input_file", position: 1 }, doc: "raw miRanda output containing target predictions" }
  output_parsed_file: { type: File, inputBinding: { prefix: "--output_file", position: 1 }, doc: "parsed and filtered results as CSV" }
  ram: { type: 'int?', default: 8, doc: "GB of RAM to allocate to the task." }
  cores: { type: 'int?', default: 2, doc: "Minimum reserved number of CPU cores for the task." }
outputs:
  array_dirs:
    type: 'File'
    outputBinding:
      glob: $(inputs.output_parsed_file)
