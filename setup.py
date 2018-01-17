from setuptools import setup

import plasmiduncover


VERSION = plasmiduncover.__version__

setup(
    name="plasmiduncover",
    version="1.0.4",
    packages=["plasmiduncover", "plasmiduncover.tools"],
    install_requires=[
        "plotly>=2.0.9",
        "biopython>=1.69",
        "termcolor>=1.1.0",
    ],
    url="https://github.com/tiagofilipe12/PlasmidCoverage",
    license="GPL3",
    author="Tiago F. Jesus",
    author_email="tiagojesus@medicina.ulisboa.pt",
    description="Script to obtain plasmid id from WGS data using bowtie2 to "
                "map",
    entry_points={
        "console_scripts": [
            "PlasmidUNCover.py = plasmiduncover.PlasmidUNCover:main",
            "diffs_json.py = plasmiduncover.tools.diffs_json:main",
            "plasmid_or_not.py = plasmiduncover.tools.plasmid_or_not:main"
        ]
    }
)
