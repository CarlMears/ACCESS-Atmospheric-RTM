import nox

# By default just run these basic lint jobs
nox.options.sessions = ["black", "flake8", "mypy"]


@nox.session
def black(session: nox.Session) -> None:
    """Check if black needs to be run"""
    session.install("black")
    session.run(
        "black", "--check", "--target-version=py38", "--diff", "access_atmosphere/"
    )


@nox.session
def flake8(session: nox.Session) -> None:
    """Run flake8"""
    session.install("flake8", "flake8-docstrings", "flake8-import-order")
    session.run("flake8", "--count", "access_atmosphere/")


@nox.session
def mypy(session: nox.Session) -> None:
    """Run mypy"""
    session.install("mypy", "numpy", "netCDF4", "cdsapi")
    session.run(
        "mypy",
        "--pretty",
        "--show-error-context",
        "access_atmosphere/",
    )


@nox.session
def mypy_full(session: nox.Session) -> None:
    """Run mypy and generate report files"""
    session.install("mypy", "numpy", "netCDF4", "cdsapi", "lxml")
    session.run(
        "mypy",
        "--pretty",
        "--show-error-context",
        "--junit-xml=mypy.junit.xml",
        "--cobertura-xml-report=.",
        "--lineprecision-report=.",
        "access_atmosphere/",
    )
