from setuptools import setup  

if __name__ == '__main__':
    with open("README.txt", "r", encoding="utf-8") as fh:
        long_description = fh.read()
    setup(name='Intercaat',
        version="3.3",
        description='This program uses a PDB file to identify the residues present in the interface between a query chain and an interacting chain(s)',
        url='https://gitlab.com/fiserlab.org/intercaat',
        author='Steve Grudman',
        author_email='steven.grudman@einsteinmed.edu',
        license='MIT',
        packages=['intercaat'],
        package_data={'intercaat': ['intercaat_config.ini'] }, 
        install_requires=[         
            'scipy',         
            'numpy',
            'pyhull'],
        long_description=long_description,
        long_description_content_type='text/markdown',
        entry_points = {
            'console_scripts': ['intercaat=intercaat.intercaat:main'],
            
        })

#python setup.py sdist
#python3 -m twine upload dist/*