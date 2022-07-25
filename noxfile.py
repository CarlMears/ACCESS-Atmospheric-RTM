import nox

nox.options.sessions = ["black", "flake8", "mypy_brief"]


@nox.session
def black(session: nox.Session) -> None:
    session.install("black")
    session.run(
        "black", "--check", "--target-version=py38", "--diff", "access_atmosphere/"
    )


@nox.session
def flake8(session: nox.Session) -> None:
    session.install("flake8", "flake8-docstrings", "flake8-import-order")
    session.run("flake8", "--count", "access_atmosphere/")


@nox.session
def mypy_brief(session: nox.Session) -> None:
    """Run mypy without extra reports"""
    session.install("mypy", "numpy")
    session.run(
        "mypy",
        "--pretty",
        "--ignore-missing-imports",
        "--show-error-context",
        "access_atmosphere/",
    )


@nox.session
def mypy_full(session: nox.Session) -> None:
    """Run mypy and generate report files"""
    session.install("mypy", "numpy", "lxml")
    session.run(
        "mypy",
        "--pretty",
        "--ignore-missing-imports",
        "--show-error-context",
        "--junit-xml=mypy.junit.xml",
        "--cobertura-xml-report=.",
        "--lineprecision-report=.",
        "access_atmosphere/",
    )
