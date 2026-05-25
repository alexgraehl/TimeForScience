#!/usr/bin/env python3

"""
### One-time setup to create a virtual environment

You'll need the python libraries for LM Studio. (LM Studio is a standalone application that you can download and run on the Mac as a regular program.)

To set up `venv` and install python's `lmstudio` library for the first time:

```sh
cd path/to/the/directory/for/image_translation_via_llm
python3 -m venv venv  # Note that the second 'venv' here is the directory name, which you could theoretically change the name of
source venv/bin/activate
pip install lmstudio
```

### Separately, you'll need to run a vision model in LM Studio.app:

Here, I'm using LM Studio. That's a standalone Mac application. You could probably also run this in Ollama, but you'd need to slightly change the API calls here.

For LM Studio, you just double click it and run it (and load a vision language model, like QwenVL or Gemma4).

### Now that LM Studio is running AND you've selected a vision language model in it

In a terminal window, activate the virtual environmen with `source`, and then run the translator:

```sh
cd path/to/the/directory/for/image_translation_via_llm
source venv/bin/activate

# Try this to get a translation into English and NOT write the output text files
python3 image_translator_calls_lmstudio.py --dry ./images/spanish/example-garfield-crop-spanish.png

# Or, say, to Japanese, for multiple files:
python3 image_translator_calls_lmstudio.py --recursive --dry --lang=Japanese ./images/spanish/ ./images/chinese/
```

"""

from collections.abc import Sequence
from pathlib import Path
from typing import Any
import argparse
import lmstudio as lms
import time

# If an image can't be processed for some reason, we'll add it here and print all problematic files at the end.
GLOBAL_FAILURES: list[str] = []

# Store the inference time for each image (image path -> inference time (sec)).
GLOBAL_TIMINGS: dict[str, float] = {}

# All recognized image file extensions. These should be lower-case!
# (The actual check is case-insensitive.) Note that the vision language model
# will also need to be able to process this format.
IMAGE_EXTENSIONS = ('.jpg', '.jpeg', '.png', '.gif', '.webp')

PROMPT_TEMPLATE = """
Transcribe the text in this image and translate it into {output_language}.
If there are multiple blocks of text, output each block separately, with a "----" divider between them.

The source language may be any language. Identify it as $LANGUAGE.

For each block of text, format your response exactly like this:
# $LANGUAGE: <original text here>
# {output_language}: <{output_language} translation here>
"""

def is_recognized_image_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in IMAGE_EXTENSIONS

def de_duplicate_items(paths: Sequence[Any]) -> Sequence[Any]:
    """De-duplicates an input and preserves order. Uses the built-in Pyton ordered-dict properties."""
    ordered_setlike: dict[Any, None] = dict.fromkeys(paths) # Relies on this being an ordered dict, by the way, which is the default these days.
    return list(ordered_setlike.keys())

def get_image_files(paths: Sequence[str], recursive: bool) -> Sequence[Path]:
    """Gathers image files from a list of paths, optionally searching directories recursively."""
    # De-duplciate the input paths, preserving order when practical
    
    image_files: list[Path] = []
    for path_str in de_duplicate_items(paths):
        path_obj = Path(path_str)
        if not path_obj.exists():
            print(f"Warning: Path does not exist: {path_str}")
            continue

        if is_recognized_image_file(path_obj):
            image_files.append(path_obj)
        elif path_obj.is_dir():
            # `glob` handles the recursion for us, so we don't need to do it ourselves. Neat!
            for file in path_obj.glob("**/*" if recursive else "*"):
                if is_recognized_image_file(file):
                    image_files.append(file)
    return de_duplicate_items(image_files)

def process_images(paths: Sequence[str], recursive: bool, dry_run: bool, overwrite: bool, output_language: str):
    """Translate text in images and optionally save the results.

    Note that the user must have already separately launched LM Studio.app and
    loaded a vision model (e.g., Gemma4).
    
    For the args, see the argparse section in `main`.
    """
    please_translate_my_image_prompt: str = PROMPT_TEMPLATE.format(output_language=output_language)

    try:
        model = lms.llm()
    except Exception as e:
        print(f"[:ERR:] Can't connect to `LM Studio.app`. Is it running (on your local machine), and did you ALSO load a model? Error: {e}")
        return

    all_img_paths: Sequence[Path] = get_image_files(paths, recursive)
    if not all_img_paths:
        print("[:ERR:] No image files were found! Double-check your input paths.")
        return

    print(f"[:FYI:] Found a total of {len(all_img_paths)} image(s). Attempting to translate them now...")

    # The suffix that we add to each translation text file. Note that this includes
    # the original image file path, including in its file type (e.g. "x.jpeg.translation.english.txt").
    translation_suffix: str = f"translation.{output_language.lower()}.txt"
    for img_path in all_img_paths:
        translated_output_path: Path = img_path.with_name(f"{img_path.name}.{translation_suffix}")  # E.g. "cool_image.translation.english.txt"
        
        # If it's a dry run, we don't care if the output file exists.
        # If the user has specified --force (overwrite), then we don't care if the file exists.
        # Otherwise, only skip the image if the translation file already exists.
        if not dry_run and (translated_output_path.exists() and not overwrite):
            print(f"[:FYI:] Skipping (since a translation exists already): `{img_path.name}`")
            continue

        print(f"[:FYI:] Translating `{img_path}`...")
        try:
            image_handle = lms.prepare_image(str(img_path.absolute()))
            chat = lms.Chat()
            
            # The message to the chat is our prompt PLUS the image.
            # Note that the prompt is (currently) the same for all images.
            chat.add_user_message(please_translate_my_image_prompt, images=[image_handle])

            # This is the slowest single step
            start_time: float = time.perf_counter()
            response = model.respond(chat)
            duration: float = time.perf_counter() - start_time
            translation: str = response.content  # We could postprocess it if we wanted, too
            GLOBAL_TIMINGS[str(img_path)] = duration

            print(f"[:FYI:] Translation of {img_path.name}:\n{translation}\n")
            if dry_run:
                print(f"[:FYI:] --dry_run is enabled, so we are not writing `{translated_output_path}`")
                continue
            
            with open(translated_output_path, "w", encoding="utf-8") as f:
                print(f"[:FYI:] Saving translation to `{translated_output_path}`...")
                f.write(translation)
            print(f"[:FYI:] (Saved translation to `{translated_output_path}`)")
            pass
        except Exception as e:
            failure_message = f"Failed to translate the file `{img_path.name}`: {e}"
            GLOBAL_FAILURES.append(failure_message)
            print(f"[:ERR:] {failure_message}")

if __name__ == "__main__":
    global_start_time: float = time.perf_counter()
    parser = argparse.ArgumentParser(description="Translates text found in images, using a vision model being run by LM Studio (which you must run separately)")
    parser.add_argument("-r", "--recursive",     action="store_true",                 help="If true, search directories recursively.")
    parser.add_argument("--dry", "--dry_run",    dest="dry_run", action="store_true", help="Do NOT write output files (but still print the translations to the terminal).")
    parser.add_argument("-f", "--force",         action="store_true",                 help="Force-overwrite any output translation files that already exist.")
    parser.add_argument("-l", "--lang",          default="English",                   help="Output language for translation (default: English).")
    
    # Accept one or more paths (files or directories) at the end of the command line.
    # Note that you just provide the paths (i.e. you do not specify "--paths").
    parser.add_argument("paths", nargs="+", help="One or more directory paths or image file paths.")
    
    args = parser.parse_args()

    process_images(args.paths, args.recursive, args.dry_run, args.force, args.lang)
    global_duration: float = time.perf_counter() - global_start_time

    total_inference_time: float = sum(GLOBAL_TIMINGS.values())
    
    if GLOBAL_FAILURES:
        print(f"[:ERR:] These {len(GLOBAL_FAILURES)} image(s) could NOT be properly translated:")
        for failure in GLOBAL_FAILURES:
            print(f"[:ERR:] * {failure}")

    # {:8.2f} is right-aligning the numbers with an 8-character padding (which should hopefully be plenty)
    print("[:OK:] Global total time  : {:8.2f}s".format(global_duration))
    print("[:OK:] Inference-only time: {:8.2f}s".format(total_inference_time))
    print("[:OK:] Non-inference time : {:8.2f}s".format(global_duration - total_inference_time))
    if len(GLOBAL_TIMINGS) > 0:
        durations = list(GLOBAL_TIMINGS.values())
        print(f"[:OK:] NUM  images handled: {len(durations):5.0f} images")
        print(f"[:OK:] NUM        problems: {len(GLOBAL_FAILURES):5.0f} problems")
        print(f"[:OK:] MEAN time per image: {sum(durations) / len(durations):8.2f}s")
        print(f"[:OK:] MAX  time per image: {max(durations):8.2f}s")
        print(f"[:OK:] MIN  time per image: {min(durations):8.2f}s")

    if not GLOBAL_FAILURES:
        print("\n[:OK:] Exiting successfully!")
    pass
