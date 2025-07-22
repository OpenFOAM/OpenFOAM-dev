from pathlib import Path
import jinja2

TEMPLATE_DIR = Path(__file__).parent / 'templates'

j2_env = jinja2.Environment(loader=jinja2.FileSystemLoader(TEMPLATE_DIR))
