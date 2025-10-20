process CONVERT_XLSX_TO_CSV {
    tag "${task.ext.prefix.id}"
    label 'process_low'

    container 'docker.io/gregorysprenger/pandas-excel:v2.2.2'

    input:
    path(spreadsheet)

    output:
    path("${task.ext.prefix.id}.csv"), emit: csv
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix
    """
    #!/usr/bin/env python3

    import platform
    import subprocess
    import pandas as pd

    batch = pd.read_excel(
        "${spreadsheet}".replace("\\\\", ""),
        engine="openpyxl",
        skiprows=1
        if pd.read_excel("${spreadsheet}", nrows=1, engine="openpyxl")
        .columns[0]
        .startswith("Run")
        else 0,
    )
    batch.to_csv("${prefix.id}.csv", index=False)

    # Output version information
    with open("versions.yml", "w") as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
    """

    stub:
    def prefix = task.ext.prefix
    """
    #!/usr/bin/env python3

    import platform
    import subprocess
    import pandas as pd

    batch = pd.read_excel(
        "${spreadsheet}".replace("\\\\", ""),
        engine="openpyxl",
        skiprows=1
        if pd.read_excel("${spreadsheet}", nrows=1, engine="openpyxl")
        .columns[0]
        .startswith("Run")
        else 0,
    )
    batch.to_csv("${prefix.id}.csv", index=False)

    # Output version information
    with open("versions.yml", "w") as f:
        f.write(f'"{subprocess.getoutput("echo ${task.process}")}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
    """
}
