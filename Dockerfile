# Use the uv base image
FROM ghcr.io/astral-sh/uv:debian

# Set working directory in the container
WORKDIR /app

# Install R and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    dirmngr \
    gnupg \
    curl \
    git \
    software-properties-common

# Add R repository and install R and pygraphviz
RUN apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev graphviz graphviz-dev\
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY . /app/

# Create and activate virtual environment using uv
RUN uv sync

# # Make sure Chainlit runs in production mode
# ENV CHAINLIT_ENV=prod

# Set environment variables for Python to prevent writing .pyc files and ensuring stdout/stderr are sent straight to terminal
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Expose the port Chainlit runs on
EXPOSE 8000

# Command to run the application
CMD ["bash","-c","cd /app && source .venv/bin/activate && chainlit run main.py --host 0.0.0.0"]