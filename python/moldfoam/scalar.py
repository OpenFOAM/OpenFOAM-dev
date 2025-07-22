import os
import re
import numpy as np
from pathlib import Path
from collections import OrderedDict
from jinja2 import Template

class OpenFOAMFieldReader:
    """Class for reading OpenFOAM field files."""
    
    def __init__(self, case_dir=None):
        """Initialize the reader with an optional case directory."""
        self.case_dir = Path(case_dir) if case_dir else None
        
    def read_field(self, field_path):
        """
        Read an OpenFOAM field file and return the header info and data.
        
        Args:
            field_path: Path to the field file
            
        Returns:
            dict: A dictionary containing header information and field data
        """
        if self.case_dir and not os.path.isabs(field_path):
            field_path = self.case_dir / field_path
            
        with open(field_path, 'r') as f:
            content = f.read()
        
        # Parse the header
        header = self._parse_header(content)
        
        # Parse the internal field
        internal_field = self._parse_internal_field(content)
        
        # Parse boundary fields
        boundary_fields = self._parse_boundary_fields(content)

        # Parse additional variables defined outside main sections
        additional_vars = self._parse_additional_variables(content)
        
        return {
            'header': header,
            'internal_field': internal_field,
            'boundary_fields': boundary_fields,
            'additional_vars': additional_vars,
            'raw_content': content,
        }
    
    def _parse_header(self, content):
        """Parse the header section of an OpenFOAM field file."""
        header = {}
        
        # Extract FoamFile dictionary
        foam_file_match = re.search(r'FoamFile\s*\{(.*?)\}', content, re.DOTALL)
        if foam_file_match:
            foam_file_content = foam_file_match.group(1)
            
            # Extract key-value pairs
            for match in re.finditer(r'\s*(\w+)\s*(.*?);', foam_file_content):
                key, value = match.groups()
                header[key] = value.strip()
        
        # Extract dimensions
        dimensions_match = re.search(r'dimensions\s*\[(.*?)\]', content)
        if dimensions_match:
            dimensions_str = dimensions_match.group(1)
            header['dimensions'] = [float(d) if '.' in d else int(d) 
                                    for d in dimensions_str.split()]
        
        return header
    
    def _parse_internal_field(self, content):
        """Parse the internal field section of an OpenFOAM field file."""
        # Check if uniform or nonuniform
        uniform_match = re.search(r'internalField\s+uniform\s+(.*?);', content)
        if uniform_match:
            value_str = uniform_match.group(1)
            # Handle scalar, vector, tensor values
            try:
                # Try parsing as a number
                return float(value_str)
            except ValueError:
                # Try parsing as a vector/tensor (remove parentheses and split)
                if '(' in value_str and ')' in value_str:
                    vector_str = value_str.strip('()')
                    return np.array([float(v) for v in vector_str.split()])
                return value_str
        
        # Find nonuniform field using a more robust approach
        internal_field_match = re.search(r'internalField\s+nonuniform\s+List<\w+>\s*\n?\s*(\d+)\s*\n?\s*\(', content)
        if not internal_field_match:
            return None
            
        size = int(internal_field_match.group(1))
        list_start = internal_field_match.end()
        
        # Find the closing parenthesis and semicolon by counting parentheses
        # This is more robust than using a simple regex with .*?
        level = 1  # Start with level 1 since we already found the opening parenthesis
        pos = list_start
        
        while level > 0 and pos < len(content):
            if content[pos] == '(':
                level += 1
            elif content[pos] == ')':
                level -= 1
            pos += 1
            
        if level != 0 or pos >= len(content):
            # Unbalanced parentheses
            print("Warning: Unbalanced parentheses in internalField section")
            return None
            
        # Check if there's a semicolon after the closing parenthesis
        semicolon_pos = content.find(';', pos)
        if semicolon_pos == -1 or not content[pos:semicolon_pos].strip() == '':
            print("Warning: No semicolon found after internalField list closing parenthesis")
            return None
            
        # Extract the data string
        data_str = content[list_start:pos-1]  # Exclude the final closing parenthesis
            
        # Remove comments
        data_str = re.sub(r'//.*?\n', '\n', data_str)
        
        # Parse data values
        values = []
        for line in data_str.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
                
            # Handle scalar values
            if '(' not in line:
                values.extend([float(v) for v in line.split()])
            else:
                # Handle vector/tensor values
                for vector_match in re.finditer(r'\((.*?)\)', line):
                    vector_str = vector_match.group(1)
                    values.append(np.array([float(v) for v in vector_str.split()]))
        
        # Verify we have the expected number of values
        if len(values) != size:
            print(f"Warning: Expected {size} values in internalField, but got {len(values)}")
        
        return np.array(values)
    
    def _parse_boundary_fields(self, content):
        """Parse the boundary fields section of an OpenFOAM field file."""
        boundary_fields = OrderedDict()
        
        # First find the boundaryField section
        boundary_section_match = re.search(r'boundaryField\s*\{', content)
        if not boundary_section_match:
            return boundary_fields
            
        # Get the start position of the boundaryField dictionary
        start_pos = boundary_section_match.end()
        
        # Extract the entire boundaryField content by properly handling nested braces
        # This is more reliable than using a simple regex
        level = 1  # Start with level 1 since we already found the opening brace
        pos = start_pos
        
        # Find the matching closing brace by counting opening and closing braces
        while level > 0 and pos < len(content):
            if content[pos] == '{':
                level += 1
            elif content[pos] == '}':
                level -= 1
            pos += 1
            
        if level != 0:
            # Unbalanced braces, return empty dictionary
            print("Warning: Unbalanced braces in boundaryField section")
            return boundary_fields
            
        # Extract the content of the boundaryField dictionary
        boundary_content = content[start_pos:pos-1]  # Exclude the final closing brace
        
        # Now parse each patch entry
        # Find all patch entries using the same approach
        pattern = re.compile(r'\s*([^\s{]+)\s*\{', re.DOTALL)
        
        # Find all patch entry starting positions
        patch_starts = [(m.group(1), m.end()) for m in pattern.finditer(boundary_content)]

        for i, (patch_name, patch_start) in enumerate(patch_starts):
            # Find the end of this patch entry by counting braces
            level = 1  # Start with level 1 for the opening brace we already found
            pos = patch_start
            
            while level > 0 and pos < len(boundary_content):
                if boundary_content[pos] == '{':
                    level += 1
                elif boundary_content[pos] == '}':
                    level -= 1
                pos += 1
                
            if level != 0:
                print(f"Warning: Unbalanced braces in patch {patch_name}")
                continue
                
            # Extract patch content (excluding braces)
            patch_content = boundary_content[patch_start:pos-1]
            
            patch_data = {}
            
            # Extract type
            type_match = re.search(r'type\s+(.*?);', patch_content)
            if type_match:
                patch_data['type'] = type_match.group(1).strip()
            
            # Extract value if present
            value_match = re.search(r'value\s+(.*?);', patch_content, re.DOTALL)
            if value_match:
                value_str = value_match.group(1).strip()
                
                # Handle uniform value
                if value_str.startswith('uniform'):
                    uniform_value = value_str.replace('uniform', '').strip()
                    
                    # Handle scalar
                    try:
                        patch_data['value'] = float(uniform_value)
                    except ValueError:
                        # Handle vector/tensor
                        if '(' in uniform_value and ')' in uniform_value:
                            vector_str = uniform_value.strip('()')
                            patch_data['value'] = np.array([float(v) for v in vector_str.split()])
                        else:
                            patch_data['value'] = uniform_value
                
                # Handle nonuniform value (List)
                elif 'nonuniform' in value_str:
                    # Find the List size and opening parenthesis
                    size_match = re.search(r'nonuniform\s+List<\w+>\s*\n?\s*(\d+)\s*\n?\s*\(', value_str)
                    if size_match:
                        size = int(size_match.group(1))
                        list_start = size_match.end()
                        
                        # Find the matching closing parenthesis
                        list_end = value_str.rfind(')')
                        if list_end == -1:
                            print(f"Warning: No closing parenthesis found for List in patch {patch_name}")
                            continue
                            
                        data_str = value_str[list_start:list_end]
                        
                        # Remove comments
                        data_str = re.sub(r'//.*?\n', '\n', data_str)
                        
                        # Parse data values
                        values = []
                        for line in data_str.strip().split('\n'):
                            line = line.strip()
                            if not line:
                                continue
                                
                            # Handle scalar values
                            if '(' not in line:
                                values.extend([float(v) for v in line.split()])
                            else:
                                # Handle vector/tensor values
                                for vector_match in re.finditer(r'\((.*?)\)', line):
                                    vector_str = vector_match.group(1)
                                    values.append(np.array([float(v) for v in vector_str.split()]))
                        
                        patch_data['value'] = np.array(values)
            
            # Extract other parameters
            for param_match in re.finditer(r'(\w+)\s+(.*?);', patch_content):
                key, value = param_match.groups()
                if key not in ['type', 'value']:
                    patch_data[key] = value.strip()
            
            boundary_fields[patch_name] = patch_data
        
        return boundary_fields
    
    def _parse_additional_variables(self, content):
        """Parse additional global variables defined outside main sections."""
        additional_vars = OrderedDict()
        
        # Create a copy of content to work with
        working_content = content
        
        # Remove main sections to avoid parsing their content as variables
        
        # Remove FoamFile section
        foam_file_match = re.search(r'FoamFile\s*\{', working_content)
        if foam_file_match:
            start_pos = foam_file_match.start()
            pos = foam_file_match.end()
            level = 1
            while level > 0 and pos < len(working_content):
                if working_content[pos] == '{':
                    level += 1
                elif working_content[pos] == '}':
                    level -= 1
                pos += 1
            working_content = working_content[:start_pos] + working_content[pos:]
        
        # Remove dimensions line
        working_content = re.sub(r'dimensions\s*\[[^\]]*\]\s*;', '', working_content)
        
        # Remove internalField section
        internal_match = re.search(r'internalField\s+', working_content)
        if internal_match:
            start_pos = internal_match.start()
            pos = internal_match.end()
            
            # Handle uniform case
            uniform_match = re.match(r'uniform\s+[^;]*;', working_content[pos:])
            if uniform_match:
                end_pos = pos + uniform_match.end()
                working_content = working_content[:start_pos] + working_content[end_pos:]
            else:
                # Handle nonuniform case - look for semicolon after balanced parentheses
                paren_level = 0
                brace_level = 0
                while pos < len(working_content):
                    char = working_content[pos]
                    if char == '(':
                        paren_level += 1
                    elif char == ')':
                        paren_level -= 1
                    elif char == '{':
                        brace_level += 1
                    elif char == '}':
                        brace_level -= 1
                    elif char == ';' and paren_level == 0 and brace_level == 0:
                        pos += 1
                        break
                    pos += 1
                working_content = working_content[:start_pos] + working_content[pos:]
        
        # Remove boundaryField section
        boundary_match = re.search(r'boundaryField\s*\{', working_content)
        if boundary_match:
            start_pos = boundary_match.start()
            pos = boundary_match.end()
            level = 1
            while level > 0 and pos < len(working_content):
                if working_content[pos] == '{':
                    level += 1
                elif working_content[pos] == '}':
                    level -= 1
                pos += 1
            working_content = working_content[:start_pos] + working_content[pos:]
        
        # Remove C++ style comments
        working_content = re.sub(r'//.*?\n', '\n', working_content)
        # Remove C style comments
        working_content = re.sub(r'/\*.*?\*/', '', working_content, flags=re.DOTALL)
        
        # Now find variable definitions in remaining content
        # Look for patterns like: variableName content;
        # We need to handle nested structures properly for the content part
        
        lines = working_content.split('\n')
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines
            if not line:
                i += 1
                continue
            
            # Look for variable definition start (word at beginning of line)
            var_match = re.match(r'^(\w+)\s+(.*)', line)
            if var_match:
                var_name = var_match.group(1)
                remaining_content = var_match.group(2)
                
                # Collect the complete variable definition until we find the terminating semicolon
                # Need to handle nested parentheses and braces
                var_lines = [line]
                paren_level = remaining_content.count('(') - remaining_content.count(')')
                brace_level = remaining_content.count('{') - remaining_content.count('}')
                
                # Check if this line already ends the definition
                if remaining_content.rstrip().endswith(';') and paren_level == 0 and brace_level == 0:
                    # Complete definition on single line
                    full_definition = line
                else:
                    # Multi-line definition - keep collecting until balanced and ends with semicolon
                    i += 1
                    while i < len(lines):
                        current_line = lines[i]
                        var_lines.append(current_line)
                        
                        # Update nesting levels
                        paren_level += current_line.count('(') - current_line.count(')')
                        brace_level += current_line.count('{') - current_line.count('}')
                        
                        # Check if we've reached the end
                        if (current_line.rstrip().endswith(';') and 
                            paren_level == 0 and brace_level == 0):
                            break
                        i += 1
                    
                    full_definition = '\n'.join(var_lines)
                
                # Store the complete variable definition
                additional_vars[var_name] = full_definition
            
            i += 1
        
        return additional_vars


class OpenFOAMFieldWriter:
    """Class for writing OpenFOAM field files."""
    
    def __init__(self, case_dir=None):
        """Initialize the writer with an optional case directory."""
        self.case_dir = Path(case_dir) if case_dir else None
        
        # Template for scalar field file
        self.scalar_template = Template('''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  {{version}}                            |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     {{foam_version}};
    format      {{format}};
    class       {{class_}};
    location    "{{location}}";
    object      {{object}};
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [{{ dimensions|join(' ') }}];
                                        
{% for name, var_content in additional_vars.items() %}
{{ var_content }}
                                        
{% endfor %}

internalField   {% if is_uniform %}uniform {{ value }}{% else %}nonuniform List<scalar>
{{ values|length }}
(
{% for val in values %}
    {{ val }}{% endfor %}
){% endif %};

boundaryField
{
{% for name, patch in boundary_fields.items() %}
    {{ name }}
    {
        type            {{ patch.type }};
        {% if 'value' in patch %}value           {% if patch.is_uniform %}uniform {{ patch.value }}{% else %}nonuniform List<scalar>
        {{ patch.values|length }}
        (
        {% for val in patch.values %}
            {{ val }}{% endfor %}
        ){% endif %};
        {% endif %}
        {% for key, val in patch.items() %}{% if key not in ['type', 'value', 'is_uniform', 'values'] %}
        {{ key }}            {{ val }};{% endif %}{% endfor %}
    }
{% endfor %}
}

// ************************************************************************* //
''')

    def write_field(self, field_data, output_path):
        """
        Write field data to an OpenFOAM field file.
        
        Args:
            field_data: Dictionary containing header info and field data
            output_path: Path where the field file will be written
        """
        if self.case_dir and not os.path.isabs(output_path):
            output_path = self.case_dir / output_path

        # Prepare header information
        header = field_data.get('header', {})
        internal_field = field_data.get('internal_field')
        boundary_fields = field_data.get('boundary_fields', {})
        additional_vars = field_data.get('additional_vars', {})
        
        # Format dimensions if present
        dimensions = header.get('dimensions', [0, 0, 0, 0, 0, 0, 0])
        
        # Check if internal field is uniform
        is_uniform = not isinstance(internal_field, np.ndarray)
        
        # Format internal field value
        if is_uniform:
            if isinstance(internal_field, np.ndarray):  # Vector/tensor
                value = f"({' '.join(str(v) for v in internal_field)})"
            else:  # Scalar
                value = str(internal_field)
            values = []
        else:
            value = None
            values = internal_field
        
        # Format boundary fields
        formatted_boundary_fields = {}
        for patch_name, patch_data in boundary_fields.items():
            formatted_patch = patch_data.copy()
            
            if 'value' in patch_data:
                patch_value = patch_data['value']
                patch_is_uniform = not isinstance(patch_value, np.ndarray) or patch_value.ndim == 1
                
                if patch_is_uniform:
                    if isinstance(patch_value, np.ndarray):  # Vector/tensor
                        formatted_patch['value'] = f"({' '.join(str(v) for v in patch_value)})"
                    else:  # Scalar
                        formatted_patch['value'] = str(patch_value)
                    formatted_patch['is_uniform'] = True
                else:
                    formatted_patch['values'] = patch_value
                    formatted_patch['is_uniform'] = False
            
            formatted_boundary_fields[patch_name] = formatted_patch
        
        # Render template
        output_content = self.scalar_template.render(
            version=header.get('version', ''),
            foam_version=header.get('version', '2.0'),
            format=header.get('format', 'ascii'),
            class_=header.get('class', 'volScalarField'),
            location=header.get('location', '0'),
            object=header.get('object', 'field'),
            dimensions=dimensions,
            additional_vars=additional_vars,
            is_uniform=is_uniform,
            value=value,
            values=values,
            boundary_fields=formatted_boundary_fields
        )
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Write file
        with open(output_path, 'w') as f:
            f.write(output_content)
        
        return output_path


class ScalarFieldProcessor:
    """High-level class for processing OpenFOAM field files."""
    
    def __init__(self, case_dir=None):
        """Initialize the processor with an optional case directory."""
        self.case_dir = Path(case_dir) if case_dir else None
        self.reader = OpenFOAMFieldReader(case_dir)
        self.writer = OpenFOAMFieldWriter(case_dir)
    
    def read_field(self, field_path):
        """Read a field file."""
        return self.reader.read_field(field_path)
    
    def write_field(self, field_data, output_path):
        """Write field data to a file."""
        return self.writer.write_field(field_data, output_path)
    
    def apply_function(self, field_data, func):
        """
        Apply a function to the internal field data.
        
        Args:
            field_data: Field data dictionary
            func: Function to apply to the internal field
            
        Returns:
            dict: Updated field data dictionary
        """
        updated_data = field_data.copy()
        internal_field = field_data['internal_field']
        
        if isinstance(internal_field, np.ndarray):
            updated_data['internal_field'] = func(internal_field)
        else:
            updated_data['internal_field'] = func(internal_field)
            
        return updated_data
    
    def combine_fields(self, field_a, field_b, operation='add'):
        """
        Combine two field datasets with a specified operation.
        
        Args:
            field_a: First field data dictionary
            field_b: Second field data dictionary
            operation: Operation to perform ('add', 'subtract', 'multiply', 'divide')
            
        Returns:
            dict: Combined field data dictionary
        """
        # Ensure both fields have the same format
        if field_a['header'].get('class') != field_b['header'].get('class'):
            raise ValueError("Field classes don't match")
            
        result = field_a.copy()
        
        # Combine internal fields
        a_internal = field_a['internal_field']
        b_internal = field_b['internal_field']
        
        if operation == 'add':
            result['internal_field'] = a_internal + b_internal
        elif operation == 'subtract':
            result['internal_field'] = a_internal - b_internal
        elif operation == 'multiply':
            result['internal_field'] = a_internal * b_internal
        elif operation == 'divide':
            result['internal_field'] = a_internal / b_internal
        else:
            raise ValueError(f"Unsupported operation: {operation}")
            
        return result
    
    def scale_field(self, field_data, factor):
        """
        Scale a field by a constant factor.
        
        Args:
            field_data: Field data dictionary
            factor: Scaling factor
            
        Returns:
            dict: Scaled field data dictionary
        """
        return self.apply_function(field_data, lambda x: x * factor)
    
    def extract_values_at_locations(self, field_data, locations):
        """
        Extract field values at specific coordinates (requires cell mapping).
        This is a placeholder function - implementation would require a mesh reader.
        
        Args:
            field_data: Field data dictionary
            locations: List of (x, y, z) coordinates
            
        Returns:
            dict: Dictionary mapping location to field value
        """
        # This would require reading the mesh to map coordinates to cell indices
        # Placeholder implementation
        raise NotImplementedError("Requires mesh reader to map coordinates to cells")
    
    def compute_statistics(self, field_data):
        """
        Compute basic statistics for a field.
        
        Args:
            field_data: Field data dictionary
            
        Returns:
            dict: Dictionary of statistics
        """
        internal_field = field_data['internal_field']
        
        if not isinstance(internal_field, np.ndarray):
            # For uniform fields, return a single value
            return {
                'min': internal_field,
                'max': internal_field,
                'mean': internal_field,
                'std': 0.0
            }
            
        return {
            'min': np.min(internal_field),
            'max': np.max(internal_field),
            'mean': np.mean(internal_field),
            'std': np.std(internal_field)
        }
