from setuptools import setup  




if __name__ == '__main__':
    with open("README.txt", "r", encoding="utf-8") as fh:
        long_description = fh.read()
    setup(name='Intercaat',
        version="1.13",
        description='This program uses a PDB file to identify the residues present in the interface between a query chain and an interacting chain(s)',
        url='https://gitlab.com/fiserlab.org/intercaat',
        author='Steve Grudman',
        author_email='steven.grudman@einsteinmed.edu',
        license='MIT',
        packages=['Intercaat'],
        package_data={'Intercaat': ['qhull/bin/*'],'Intercaat': ['intercaat_config.ini'] }, 
        install_requires=[         
            'scipy',         
            'numpy',
            'pyhull'],
        long_description=long_description,
        long_description_content_type='text/markdown',
        entry_points = {
            'console_scripts': ['intercaat=Intercaat.intercaat:main',"intercaat_setup=Intercaat.intercaat:setup" ],
            
        })

#command = "curl www.qhull.org/download/qhull-2020-src-8.0.2.tgz --output qhull-2020-src-8.0.2.tgz && tar -xvf qhull-2020-src-8.0.2.tgz && cd qhull-2020.2 &&  make &&  export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH &&  make install && cp bin/qvoronoi ../Intercaat/qhull/bin && cd .. && rm -r qhull-2020.2"
#process = subprocess.run(command, shell=True)
#python setup.py sdist
#python3 -m twine upload dist/*