import asyncio
from mcp_agent.core.fastagent import FastAgent

fast = FastAgent("FastAgent Example")


@fast.agent(instruction="You are agent who helps with analyzing csv files.", servers=['FileSystem', 'SamiAI'])
async def main():
    async with fast.run() as agent:
        # response = await agent("Can you list all the files available in the /home/niklaus/projects/biochem/sami_ai/client and dump it in a markdown file?")
        # print(response)
        await agent()


# @fast.agent(
#     instruction="You are an agent who helps with analyzing CSV files and running normalization workflows.",
#     servers=["FileSystem", "SamiAI"],
# )
# async def main():
#     async with fast.run() as agent:
#         # Example of how to run the agent with a specific prompt
#         await agent(
#             "I need to normalize the brain_metabolomics.csv file and visualize results for the TG(16:0_18:1_18:2) molecule. Can you help me with that?"
#         )


if __name__ == "__main__":
    asyncio.run(main())
