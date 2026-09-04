"""
The function of this script is to render the Jinja2 templates in the current
directory into input files for sphinx. It introspects the OpenMMForceFields
Python module to find all of the classes and formats them for inclusion into the
templates.
"""

from os.path import dirname, join, splitext, basename
from glob import glob
import inspect

import jinja2
import openmmforcefields.generators


def fullname(klass):
    return klass.__module__ + "." + klass.__name__


def template_variables():
    """Create the data structure available to the Jinja2 renderer when
    filling in the templates.
    """
    data = {
        "functions": [],
        "classes": [],
    }
    for name, obj in inspect.getmembers(openmmforcefields.generators, predicate=inspect.isfunction):
        data["functions"].append("openmmforcefields.generators." + name)
    for name, obj in inspect.getmembers(openmmforcefields.generators, predicate=inspect.isclass):
        data["classes"].append("openmmforcefields.generators." + name)
    return data


def main():
    here = dirname(__file__)
    templateLoader = jinja2.FileSystemLoader(here)
    templateEnv = jinja2.Environment(loader=templateLoader)
    data = template_variables()

    for template_fn in map(basename, glob(join(here, "*.jinja2"))):
        output_fn = splitext(template_fn)[0]
        print(f"Rendering {template_fn} to {output_fn}...")
        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, "w") as f:
            f.write(output_text)


if __name__ == "__main__":
    main()
