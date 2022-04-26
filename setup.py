from setuptools import setup

setup(name='funniest',
      version='0.1',
      description='The funniest joke in the world',
      url='',
      author='Steve Grudman',
      author_email='steven.grudman@einsteinmed.edu',
      license='MIT',
      packages=['Intercaat'],
      zip_safe=False, 
      entry_points = {
        'console_scripts': ['intercaat=Intercaat.intercaat:main'],
    })


