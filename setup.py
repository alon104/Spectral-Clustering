from setuptools import setup, Extension


setup(
	name='mykmeanssp',
	version='1.0',
	description='a spectral kmeans algorithm',
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['spkmeansmodule.c','spkmeans.c'],
        )
    ]
)
