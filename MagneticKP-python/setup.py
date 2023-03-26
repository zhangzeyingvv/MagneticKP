from setuptools import setup, find_packages

VERSION = '1.0.0'
DESCRIPTION = 'MagneticKP'
LONG_DESCRIPTION = 'A package that construct k.p model'

setup(
    name="magnetickp",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author="Zhang Zeying",
    author_email="zhangzeyingvv@gmail.com",
    license='GNU',
    packages=find_packages(),
    package_data={"": ["*.pkl"]},
    #install_requires=['sympy>=1.11.1','itertools','os','pickle'],
    install_requires=['sympy>=1.11.1'],
    keywords='kp model',
    classifiers= [
        'License :: GNU',
        "Programming Language :: Python :: 3",
    ]
)

