from tools.integration import IntegrationInput, integration_tool
import asyncio

async def main():
    input = IntegrationInput(
        file_path1="/home/niklaus/projects/biochem/sami_ai/clustering_052825_224640.h5ad",
        file_path2="/home/niklaus/projects/biochem/sami_ai/clustering_052825_225100.h5ad"
    )

    results = await integration_tool(input, session_id="sessio_id")
    print(results)

if __name__ == "__main__":
    asyncio.run(main())