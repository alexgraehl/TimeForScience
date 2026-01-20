#!/usr/bin/env python3

"""
Notes for using Ollama as a photo/screenshot renamer

* Run `ollama serve` in one terminal window
* In another one, run `ollama list` to list all models, then try like `ollama mixtral` or whatever to actually run the model. You'll then see the server terminal window display a bunch of text zooming by with a bunch of information as the server loads the model.

* Now you have:

  1) SERVER terminal window (no need to interact with this)
  2) CLIENT terminal window

We need to install a vision model. Let's get Qwen 2.5:
* Start `ollama serve` (on the command line). This must be running the whole time (in its own terminal window), you should never need to quit it!
* Then, in ANOTHER window, run this command.
   > ollama pull qwen2.5vl:7b
* (7B parameter model fits on an 'average' 2026 Mac). Looks like it's 6 GB.)
* We'll also need ollama for python, so `pip3 install ollama`

To set up `venv` and install python's `ollama` library for the first time:

```sh
python3 -m venv venv  # Note that the second 'venv' is the directory name
source venv/bin/activate
pip install ollama
```

And then to actually run this script, just run:

```sh
source venv/bin/activate
python3 screenshot_classifier_ollama.py
```

"""

import argparse
import os
import re
import textwrap
import time

import ollama

KNOWN_FILE_EXTENSIONS: tuple[str, ...] = ('.png', '.jpg', '.jpeg', '.webp')

# --- CONFIGURATION ---
DEFAULT_IMAGE_DIR="./images"
MODEL = "qwen2.5vl:7b"  # Run "ollama list" to get a list of models that we know about.
PROMPT = "Describe this image in 3-10 words suitable for a filename. Only provide a description, no punctuation."

"""Very basic sanitization.
Spaces become underscores. Hyphens are allowed. Most non-obviously-alphanumeric characters become an underscore.
A couple of special behaviors: changes "%" to "percent".
Doesn't change the case. Trims leading/trailing whitespace.
So e.g. "  Snake5%ç ç " becomes "Snake5Percent__"
"""
def sanitize_filename(input_filename: str):
    input_filename = input_filename.strip().replace("%", "percent").replace("$", "dollar")
    input_filename = input_filename.strip().replace("%", "percent")
    # Replace anything else "weird" with a suitable simple character. Don't care about upper/lower case.
    return re.sub(r"[^-A-Za-z0-9._]", "_", input_filename)

def rename_images(img_dir: str):
    if not os.path.exists(img_dir):
        raise IOError(f"Can't find the input image directory `{img_dir}`! Maybe it doesn't exist?")

    n_with_name_suggestion: int = 0
    n_failed: int = 0
    n_skipped: int = 0
    for filename in os.listdir(img_dir):
        if not filename.lower().endswith(KNOWN_FILE_EXTENSIONS):
            print(f"Skipping this file of a type that we don't know what to do with: {filename}")
            continue # Skip it!

        image_path = os.path.join(img_dir, filename)
        
        print(f"Processing {filename}... (note: NOT reading recursively yet)")
        
        try:
            llm_start_time = time.perf_counter()
            response = ollama.chat(
                model=MODEL,
                messages=[{
                    'role': 'user',
                    'content': PROMPT,
                    'images': [image_path]
                }]
            )
            llm_end_time = time.perf_counter()
            
            original_description = response['message']['content']
            clean_desc = sanitize_filename(original_description)
            
            file_extension = os.path.splitext(filename)[1].lower()  # e.g. ".jpg" or whatever. Includes the dot, notably.
            new_filename = f"{clean_desc}{file_extension}"
            new_path: str = os.path.join(img_dir, new_filename)
            
            # If there's a name collision, try incrementing a counter. Blow up if we need more than 1000 of these. Not robust!
            counter = 1
            while os.path.exists(new_path):
                new_path = os.path.join(img_dir, f"{clean_desc}_{counter}{file_extension}")
                counter += 1
                if counter > 1000:
                    raise IOError(f"How are we STILL getting name collisions? Check {new_path}.")

            # os.rename(image_path, new_path) # Dangerous!
            print(textwrap.dedent(
                f"""
                ## Ollama thought for {llm_end_time - llm_start_time:.1f} seconds, and wants to perform this renaming:
                {filename} --> {os.path.basename(new_path)}
                ## Ollama originally had this to say: {original_description}
                """))
            n_with_name_suggestion += 1

        except Exception as e:
            print(textwrap.dedent(
                f"""
                Error processing {filename}: {e}.
                Is ollama running already? If not, try running `ollama serve` and make sure you have the expected model (`${MODEL}`)
                """))
            n_failed += 1

    return n_with_name_suggestion, n_failed, n_skipped

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A renamer.")
    parser.add_argument("--in_dir", type=str, default="", help="The input directory. Will be read, but NOT recursively.")
    args = parser.parse_args()

    input_dir = DEFAULT_IMAGE_DIR if (args.in_dir == "") else args.in_dir

    overall_start_time = time.perf_counter()
    n_good, n_bad, n_skipped = rename_images(input_dir)
    overall_end_time = time.perf_counter()

    elapsed = overall_end_time - overall_start_time
    print(f"""\
          Suggested renaming options for {n_good} file(s) in {elapsed:.1f} seconds.
          Found {n_skipped} file(s) that didn't seem like valid image files.
          (Failed to rename {n_bad} file(s)).
          """)

