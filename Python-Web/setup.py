from setuptools import setup, find_packages

VERSION = '1.2.0'
DESCRIPTION = 'MagneticKP'
LONG_DESCRIPTION = 'A package that construct k.p model'
MSGCOREP_DEPENDENCY = 'msgcorep>=0.1.0a0,<0.2'
FASTAPI_DEPENDENCY = 'fastapi>=0.100,<1'
NUMPY_DEPENDENCY = 'numpy>=1.23,<3'
MATPLOTLIB_DEPENDENCY = 'matplotlib>=3.7,<4'
SCIPY_DEPENDENCY = 'scipy>=1.9,<2'

setup(
    name="magnetickp",
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    author="Zhang Zeying",
    author_email="zhangzeyingvv@gmail.com",
    url="https://github.com/zhangzeyingvv/MagneticKP",
    license="GPL-3.0-or-later",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "magnetickp.web": [
            "static/*",
            "static/help/*",
            "static/help/assets/*",
            "static/vendor/katex/*",
            "static/vendor/katex/fonts/*",
        ]
    },
    install_requires=[
        'sympy>=1.11.1',
        NUMPY_DEPENDENCY,
        SCIPY_DEPENDENCY,
        MATPLOTLIB_DEPENDENCY,
        FASTAPI_DEPENDENCY,
        "uvicorn[standard]>=0.23,<1",
        MSGCOREP_DEPENDENCY,
    ],
    extras_require={
        "test": [
            "pytest>=7,<9",
            "httpx>=0.24,<1",
            "gmpy2>=2.1,<3",
        ],
    },
    entry_points={
        "console_scripts": [
            "magnetickp-web=magnetickp.web.__main__:main",
        ]
    },
    keywords='kp model',
    python_requires=">=3.9",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
    ]
)
