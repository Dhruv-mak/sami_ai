# SAMI AI - Biochemical Data Analysis Assistant 🧪🔬🤖

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)


SAMI AI is an intelligent assistant designed to help mass spectrometry experts and biochemical researchers analyze complex datasets without writing code. Through a natural language chat interface, users can perform data analysis, visualization, and interpretation.

## Features

- **Interactive Chat Interface**: Communicate with SAMI AI using simple English commands
- **Data Processing**: Load and normalize CSV data files from various omics experiments
- **Visualization**: Generate plots and visual representations of molecules across samples
- **Marker Analysis**: Find significant markers that differentiate between clusters or regions
- **Clustering**: Perform clustering and other statistical analyses on biochemical data
- **Pipeline Creation**: Build and execute analysis workflows for repeated tasks
- **Multi-omics Support**: Work with metabolomics, lipidomics, glycomics and other data types

## Quick Start

### Using Docker (Recommended)

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/sami_ai.git
   cd sami_ai
   ```

2. Build and run with Docker:
   ```
   docker build -t sami-ai .
   docker run -p 8000:8000 -it sami-ai
   ```

3. Open your browser and navigate to http://localhost:8000

### Local Installation

1. Ensure you have Python 3.13+ and [UV](https://github.com/astral-sh/uv) installed

2. Clone and set up the repository:
   ```
   git clone https://github.com/yourusername/sami_ai.git
   cd sami_ai
   uv sync
   source .venv/bin/activate
   ```

3. Start the Chainlit server:
   ```
   chainlit run main.py
   ```

## Project Structure

```
sami_ai/
├── SAMI/                # Core analysis modules
├── tools/               # Backend tools for data processing
├── client/              # Client-side components and sample data
├── data/                # Dataset storage
├── .chainlit/           # Chainlit configuration
├── main.py              # Application entry point
├── Dockerfile           # Container definition
└── pyproject.toml       # Project dependencies
```

## Usage Examples

- "Load my brain_metabolomics.csv file and normalize it"
- "Show me the distribution of molecule X across samples"
- "Find markers that differentiate between clusters 1 and 2"
- "Create a heatmap of the top 20 molecules in my data"
- "Help me explore the lipidomics data from region A"

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📜 License

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License - see LICENSE(./LICENSE) file for details.

It includes and adapts code originally published by [Xin Ma](https://github.com/XinBiostats) in repository [SAMI](https://github.com/XinBiostats/SAMI) under the same license.

### Changes Made:
- Modified function definitions to make adaptable for the chatbot.