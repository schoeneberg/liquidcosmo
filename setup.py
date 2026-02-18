from setuptools import setup, find_packages

with open("README.md","r") as rfile:
  readme = rfile.read()

requirements = []
setup(name = "liquidcosmo",
      version= "0.3.8",
      author = "Nils Schöneberg",
      description="Cosmological analyses made easy",
      long_description = readme,
      long_description_content_type="text/markdown",
      url="https://github.com/schoeneberg/liquidcosmo/",
      #py_modules=['liquidcosmo','src.chain','src.folder','src.matplotlib_defaults'],
      packages=['liquidcosmo'],
      install_requires=requirements,
      classifiers=[
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      ],
)
