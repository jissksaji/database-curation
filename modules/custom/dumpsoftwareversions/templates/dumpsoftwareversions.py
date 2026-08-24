#!/usr/bin/env python3

import platform
from textwrap import dedent

import yaml


def make_versions_html(versions):
    rows = []
    for process, tools in sorted(versions.items()):
        for index, (tool, version) in enumerate(sorted(tools.items())):
            rows.append(
                f"<tr><td><samp>{process if index == 0 else ''}</samp></td>"
                f"<td><samp>{tool}</samp></td><td><samp>{version}</samp></td></tr>"
            )
    return dedent("""\\
        <table class=\"table\" style=\"width:100%\">
          <thead><tr><th>Process</th><th>Software</th><th>Version</th></tr></thead>
          <tbody>{}</tbody>
        </table>""").format("\\n".join(rows))


with open("$versions") as handle:
    observed_versions = yaml.safe_load(handle) or {}

this_module = {
    "${task.process}": {
        "python": platform.python_version(),
        "pyyaml": yaml.__version__,
    }
}
observed_versions.update(this_module)

by_module = {}
for process, tools in observed_versions.items():
    by_module.setdefault(process.split(":")[-1], tools)
by_module["Workflow"] = {
    "Nextflow": "$workflow.nextflow.version",
    "$workflow.manifest.name": "$workflow.manifest.version",
}

with open("software_versions.yml", "w") as handle:
    yaml.safe_dump(by_module, handle, sort_keys=True)
with open("software_versions_mqc.yml", "w") as handle:
    yaml.safe_dump(
        {
            "id": "software_versions",
            "section_name": "${workflow.manifest.name} software versions",
            "plot_type": "html",
            "description": "Software versions collected during this run.",
            "data": make_versions_html(by_module),
        },
        handle,
        sort_keys=False,
    )
with open("versions.yml", "w") as handle:
    yaml.safe_dump(this_module, handle, sort_keys=True)
