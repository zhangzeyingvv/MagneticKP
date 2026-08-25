"""Run the local MagneticKP Web application."""

import os


def main():
    import uvicorn

    uvicorn.run(
        "magnetickp.web.app:app",
        host=os.environ.get("MAGNETICKP_WEB_HOST", "127.0.0.1"),
        port=int(os.environ.get("MAGNETICKP_WEB_PORT", "8000")),
        reload=False,
    )


if __name__ == "__main__":
    main()
