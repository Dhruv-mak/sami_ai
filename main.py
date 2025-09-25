import logging
import json
from openai import AsyncOpenAI
import chainlit as cl
from typing import Dict, Any
from chainlit.input_widget import Select, Slider, Switch
import json
import os
import base64
from db import initialize_json
from tools.types import ToolResult, ToolResultType

# from chainlit import ClientSession
from tools.normalisation import normalisation_tool, NormalizationInput
from tools.clustering import ClusteringInput, clustering_tool
from tools.integration import IntegrationInput, integration_tool
from tools.markers import MarkersInput, markers_tool
from tools.pathways import PathwayInput, pathway_tool
from tools.pathway_viz import PathwayVizInput, pathway_viz_tool

# Configure the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("app.log"),
    ],
)


def encode_image(image_path: str) -> str:
    """Encode an image to base64."""
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode("utf-8")
    return encoded_string


logger = logging.getLogger(__name__)

client = AsyncOpenAI(
    api_key=os.getenv("API_KEY"),
    base_url="https://api.ai.it.ufl.edu/v1",
)

# Instrument the OpenAI client
cl.instrument_openai()


@cl.on_chat_start
async def start():
    logger.info("Chat started")
    settings = await cl.ChatSettings(
        [
            Select(
                id="model",
                label="LLM - Model",
                values=[
                    "gpt-4.1",
                    "gpt-4o",
                    "gpt-4.1-nano",
                    "gpt-4.1-mini",
                    "claude-3.5-sonnet-v2",
                    "mistral-small-3.1",
                    "claude-3.7-sonnet",
                    "o3-mini-medium",
                    "o3-mini-high",
                    "claude-3.7-sonnet-thinking",
                    "o3",
                    "claude-3.5-sonnet",
                    "o3-mini",
                    "o1",
                ],
                initial_index=0,
            ),
            Slider(
                id="temperature",
                label="LLM - Temperature",
                min=0,
                max=2,
                step=0.1,
                initial=0.7,
            ),
            Switch(id="stream", label="Stream Tokens", initial=True),
        ]
    ).send()
    cl.user_session.set("settings", settings)
    cl.user_session.set(
        "message_history",
        [
            {
                "role": "system",
                "content": "You are a helpful AI assistant specializing in biochemical data analysis. Your job is to minimize the input from the user so try to use default values unless specified otherwise for the functions."
                "You can help users process data through pipelines including normalization, "
                "visualization, clustering, and other analysis steps. "
                "You have access to MCP tools for executing these operations.",
            }
        ],
    )

    await cl.Message(
        content="Welcome to the Biochemical Data Analysis Assistant! I can help you process and analyze your data using various tools. "
        "You can ask me to:\n"
        "1. Load and normalize CSV data files\n"
        "2. Visualize specific molecules\n"
        "3. Perform clustering or other analyses\n"
        "4. Create multi-step pipelines for complex workflows"
    ).send()

    # initialize Database for the session
    session_id = cl.user_session.get("id")
    initialize_json(session_id)

    logger.info(f"User session id: {session_id}")


@cl.step(type="tool")
async def execute_tool(tool_name: str, tool_input: Dict[str, Any]) -> list[ToolResult]:
    """Execute an MCP tool with the given input."""
    logging.info(f"Executing tool: {tool_name}")
    logging.info(f"Tool input: {tool_input}")

    try:
        if tool_name == "normalisation_tool":
            logger.info(f"Normalisation tool input: {tool_input}")
            model = NormalizationInput(**tool_input)
            return await normalisation_tool(model, cl.user_session.get("id"))

        elif tool_name == "clustering_tool":
            try:
                model = ClusteringInput(**tool_input)
                return await clustering_tool(model, cl.user_session.get("id"))

            except Exception as e:
                logger.error(f"Error in clustering tool: {str(e)}")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        error=f"Error in clustering tool: {str(e)}",
                    )
                ]
        elif tool_name == "integration_tool":
            try:
                model = IntegrationInput(**tool_input)
                return await integration_tool(model, cl.user_session.get("id"))

            except Exception as e:
                logger.error(f"Error in integration tool: {str(e)}")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        error=f"Error in integration tool: {str(e)}",
                    )
                ]
        elif tool_name == "markers_tool":
            try:
                model = MarkersInput(**tool_input)
                return await markers_tool(model, cl.user_session.get("id"))

            except Exception as e:
                logger.error(f"Error in markers tool: {str(e)}")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        error=f"Error in markers tool: {str(e)}",
                    )
                ]
        elif tool_name == "pathway_tool":
            try:
                model = PathwayInput(**tool_input)
                return await pathway_tool(model, cl.user_session.get("id"))

            except Exception as e:
                logger.error(f"Error in pathway tool: {str(e)}")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        error=f"Error in pathway tool: {str(e)}",
                    )
                ]
        elif tool_name == "pathway_viz_tool":
            try:
                model = PathwayVizInput(**tool_input)
                return await pathway_viz_tool(model, cl.user_session.get("id"))

            except Exception as e:
                logger.error(f"Error in pathway visualization tool: {str(e)}")
                return [
                    ToolResult(
                        type=ToolResultType.text,
                        error=f"Error in pathway visualization tool: {str(e)}",
                    )
                ]
        else:
            raise ValueError(f"Tool '{tool_name}' not recognized.")
    except Exception as e:
        return [
            ToolResult(
                type=ToolResultType.text,
                error=f"Error executing tool '{tool_name}': {str(e)}",
            )
        ]


@cl.on_message
async def on_message(message: cl.Message):
    message_history = cl.user_session.get("message_history", [])
    message_history.append({"role": "user", "content": message.content})

    try:
        # Initial message for the first assistant response
        initial_msg = cl.Message(content="")
        await initial_msg.send()

        chat_params = {**cl.user_session.get("settings")}
        with open("tools.json") as file:
            tools = json.load(file)
        chat_params["tools"] = tools
        chat_params["tool_choice"] = "auto"
        logging.info(f"Tools passed: {json.dumps(tools)}")

        stream = await client.chat.completions.create(
            messages=message_history, **chat_params
        )

        initial_response = ""
        tool_calls = []

        async for chunk in stream:
            delta = chunk.choices[0].delta

            if token := delta.content or "":
                initial_response += token
                await initial_msg.stream_token(token)

            if delta.tool_calls:
                for tool_call in delta.tool_calls:
                    tc_id = tool_call.index
                    if tc_id >= len(tool_calls):
                        tool_calls.append({"name": "", "arguments": ""})

                    if tool_call.function.name:
                        tool_calls[tc_id]["name"] = tool_call.function.name

                    if tool_call.function.arguments:
                        tool_calls[tc_id]["arguments"] += tool_call.function.arguments

        # First, update message history with the initial response
        if initial_response.strip():
            message_history.append({"role": "assistant", "content": initial_response})

        # Process tool calls if any
        if tool_calls:
            for tool_call in tool_calls:
                tool_name = tool_call["name"]
                try:
                    tool_args = json.loads(tool_call["arguments"])

                    # Generate a unique tool call ID
                    tool_call_id = f"call_{len(message_history)}"

                    # Add the tool call to message history
                    message_history.append(
                        {
                            "role": "assistant",
                            "content": None,
                            "tool_calls": [
                                {
                                    "id": tool_call_id,
                                    "type": "function",
                                    "function": {
                                        "name": tool_name,
                                        "arguments": tool_call["arguments"],
                                    },
                                }
                            ],
                        }
                    )

                    with cl.Step(name=f"Executing tool: {tool_name}", type="tool"):
                        tool_results = await execute_tool(tool_name, tool_args)

                    tool_result_content = []
                    for result in tool_results:
                        # Handle different result types
                        if result.type == ToolResultType.text:
                            tool_result_content.append(
                                {
                                    "type": "text",
                                    "text": result.content,
                                }
                            )
                            await cl.Message(
                                content=f"Tool Result from {tool_name}:\n{result.content}",
                                author="Tool",
                            ).send()
                        elif result.type == ToolResultType.image:
                            image_path = result.content
                            tool_result_content.append(
                                {
                                    "type": "image_url",
                                    "image_url": {
                                        "url": f"data:image/png;base64,{encode_image(image_path)}"
                                    },
                                }
                            )
                            await cl.Message(
                                content=result.desc,
                                author="Tool",
                                elements=[
                                    cl.Image(
                                        path=image_path,
                                        caption=result.desc,
                                    )
                                ],
                            ).send()
                        elif result.type == ToolResultType.pyplot:
                            # Will implement pyplot handling later
                            pass
                        elif result.type == ToolResultType.plotly:
                            # Will implement plotly handling later
                            pass
                    if all(x.type == ToolResultType.text for x in tool_results):
                        message_history.append(
                            {
                                "role": "tool",
                                "tool_call_id": tool_call_id,
                                "content": tool_result_content,
                            }
                        )
                    else:
                        text_tool_result_content = []
                        for msg in tool_result_content:
                            if msg["type"] == "text":
                                text_tool_result_content.append(msg)
                        message_history.append(
                            {
                                "role": "tool",
                                "tool_call_id": tool_call_id,
                                "content": [
                                    {
                                        "type": "text",
                                        "text": "This is the text content of result, next images are from user role but are generated by this same tool.",
                                    }
                                ]
                                + text_tool_result_content,
                            }
                        )
                        image_tool_result_content = []
                        for msg in tool_result_content:
                            if msg["type"] == "image_url":
                                image_tool_result_content.append(msg)
                        message_history.append(
                            {
                                "role": "user",
                                "content": [
                                    {
                                        "type": "text",
                                        "text": "Since Openapi's tool role doesn't allow for image_url \
                                        in the message history content I'm sending it as a user. So, assume\
                                        that this is a tool result",
                                    }
                                ]
                                + image_tool_result_content,
                            }
                        )

                    # Create a new message for the follow-up response
                    follow_up_msg = cl.Message(content="")
                    await follow_up_msg.send()

                    # Stream the follow-up response
                    follow_up_stream = await client.chat.completions.create(
                        messages=message_history, **cl.user_session.get("settings")
                    )

                    follow_up_text = ""
                    async for chunk in follow_up_stream:
                        if token := chunk.choices[0].delta.content or "":
                            follow_up_text += token
                            await follow_up_msg.stream_token(token)

                    # Add the follow-up response to message history
                    message_history.append(
                        {"role": "assistant", "content": follow_up_text}
                    )

                except Exception as e:
                    error_msg = f"Error executing tool {tool_name}: {str(e)}"
                    error_message = cl.Message(content=error_msg)
                    await error_message.send()
                    logging.error(error_msg)

        cl.user_session.set("message_history", message_history)

    except Exception as e:
        error_message = f"Error: {str(e)}"
        await cl.Message(content=error_message).send()
        logging.error(f"Error in on_message: {str(e)}")


if __name__ == "__main__":
    from chainlit.cli import run_chainlit

    run_chainlit(__file__)
