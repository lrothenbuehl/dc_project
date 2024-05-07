# Pure chatgpt code. i cant be bothered to write my own paraser because of some stupid versioning thing.
import nbformat
import sys
import warnings

def notebook_to_script(notebook_path):
    # Load the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = nbformat.read(f, as_version=4)
    
    # Normalize the notebook
    notebook = nbformat.v4.upgrade(notebook)
    
    # Initialize script content
    script_content = ''

    # Extract code cells and concatenate into script content
    for cell in notebook.cells:
        if cell.cell_type == 'code':
            script_content += ''.join(cell.source) + '\n\n'

    # Write script content to a Python file
    script_filename = notebook_path.replace('.ipynb', '.py')
    with open(script_filename, 'w', encoding='utf-8') as f:
        f.write(script_content)

    print(f'Successfully converted notebook to script: {script_filename}')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python notebook_to_script.py <notebook_file.ipynb>')
        sys.exit(1)

    notebook_path = sys.argv[1]
    
    # Suppress MissingIDFieldWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        notebook_to_script(notebook_path)
