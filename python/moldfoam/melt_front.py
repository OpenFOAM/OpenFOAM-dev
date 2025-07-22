import os
import numpy as np
from typing import NamedTuple, Any
from pathlib import Path
from moldfoam.scalar import ScalarFieldProcessor


class FieldData(NamedTuple):
    alpha: np.ndarray
    T: np.ndarray


class FieldDataTemplate(NamedTuple):
    name: str
    field: Any
    dimensions: list[int]


def extract_data(case_dir):
    """
    Extract alpha field data from all time directories in the case.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        
    Returns:
        dict: Dictionary mapping time values to alpha field data
    """
    processor = ScalarFieldProcessor(case_dir)
    
    # Get time directories sorted numerically
    time_dirs = sorted([d for d in os.listdir(case_dir) 
                     if d.replace('.', '', 1).isdigit()], 
                     key=lambda x: float(x))

    # Exclude 9999 time step if it exists (this is the post-processing time stamp)
    if '9999' in time_dirs:
        time_dirs.remove('9999')
    
    # Dictionary to store alpha values by time
    data_by_time = {}

    # Read alpha field for each time step
    print(f"Reading field data from {len(time_dirs)} time directories...")
    for time_dir in time_dirs:
    
        alpha_file = Path(time_dir) / "alpha.melt"
        T_file = Path(time_dir) / "T"
        try:
            # Read the alpha field
            alpha_field = processor.read_field(str(alpha_file))
            T_field = processor.read_field(str(T_file))
            
            # Store only the internal field data for efficiency
            data_by_time[time_dir] = FieldData(
                alpha=alpha_field['internal_field'],
                T=T_field['internal_field'],
            )
            
            if len(data_by_time) == 1:
                # Store the first alpha field completely to use as a template later
                alpha_template = FieldDataTemplate(
                    name="meltFrontTime",
                    field=alpha_field,
                    dimensions=[0, 0, 1, 0, 0, 0, 0],  # Seconds
                )
                T_template = FieldDataTemplate(
                    name="meltFrontTemp",
                    field=T_field,
                    dimensions=[0, 0, 0, 1, 0, 0, 0],  # Kelvin
                )

                data_by_time['template'] = FieldData(
                    alpha=alpha_template, T=T_template
                )
        except Exception as e:
            print(f"Error reading {alpha_file}: {str(e)}")
    
    print(f"Successfully read data from {len(data_by_time)-1} time steps")
    return data_by_time


def cross_time(alpha_curr, alpha_prev, t_curr, t_prev, threshold):
    t_curr, t_prev = float(t_curr), float(t_prev)
    # Avoid division by zero
    if (t_curr - t_prev) == 0:
        return t_curr
    alpha_slope = (alpha_curr - alpha_prev) / (t_curr - t_prev)
    delta_t = (threshold - alpha_prev) / alpha_slope
    return t_prev + delta_t


def fill_uniform(data, times, field_name):
    """Fill uniform fields with cell values."""

    # Get number of cells from the first time step that has a spatially varying field
    uniform_values = []
    for i in range(len(times)):
        field = getattr(data[times[i]], field_name)
        if np.isscalar(field):
            uniform_values.append(field)
        else:
            break

    if i == len(times):
        raise ValueError(f"No non-uniform field found for {field_name} in the provided times.")

    num_cells = len(field)

    # Backfill all the previous data with uniform values
    for i in range(len(uniform_values)):
        data[times[i]] = data[times[i]]._replace(
            **{field_name: np.full(num_cells, uniform_values[i])}
        )

    return data, num_cells



def calc_melt_front(alpha_by_time, threshold=0.5):
    """
    Calculate the time when alpha crosses the threshold for each cell.
    
    Args:
        alpha_by_time: Dictionary mapping time values to alpha field data
        
    Returns:
        numpy.ndarray: Array of melt front times for each cell
    """
    # Remove the template entry before sorting times
    data = {k: v for k, v in alpha_by_time.items() if k != 'template'}
    
    # Get sorted time steps
    times = sorted(data.keys())
    
    if not times:
        raise ValueError("No time steps found with valid alpha data")

    data, num_cells = fill_uniform(data, times, 'alpha')
    data, num_cells = fill_uniform(data, times, 'T')

    print(f"Calculating melt front data for {num_cells} cells...")

    # Initialize melt_front_time with a large value
    melt_front_time = np.ones(num_cells) * float('inf')
    melt_front_temp = data[times[0]].T.copy()

    # For each time step
    for i in range(1, len(times)):
        t_curr = times[i]
        t_prev = times[i-1]
        alpha_curr = data[t_curr].alpha
        T_curr = data[t_curr].T
    
        if i == 1:
            alpha_prev = np.zeros_like(alpha_curr)
            T_prev = data[t_curr].T  # Avoid uniform field at t=0
        else:
            alpha_prev = data[t_prev].alpha
            T_prev = data[t_prev].T
    
        # Find cells where alpha crosses threshold between these time steps
        crossings = np.where((alpha_prev < threshold) & (alpha_curr >= threshold))[0]
    
        for cell in crossings:
            # Linear interpolation to find more precise crossing time
            crossing_time = cross_time(
                alpha_curr[cell], alpha_prev[cell], t_curr, t_prev, threshold
            )
            crossing_temp = cross_time(
                alpha_curr[cell], alpha_prev[cell], T_curr[cell], T_prev[cell], threshold
            )

            # Update melt_front_time if this is the first crossing
            if melt_front_time[cell] == float('inf'):
                melt_front_time[cell] = crossing_time
                melt_front_temp[cell] = crossing_temp

    # Count cells that were filled
    filled_cells = np.sum(melt_front_time < float('inf'))
    print(f"Found melt front times for {filled_cells} out of {num_cells} cells")

    # Set melt front time to the last time step for cells that never crossed the threshold
    # unfilled = np.where(np.isinf(melt_front_time))[0]
    # melt_front_time[unfilled] = float(times[-1])
    # melt_front_temp[unfilled] = data[times[-1]].T[unfilled]

    return melt_front_time, melt_front_temp

def write_melt_front_field(case_dir, melt_front_time, template, times):
    """
    Write the melt front time field to an OpenFOAM field file.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        melt_front_time: Array of melt front times
        template: Template fields to use for metadata
        times: List of time steps used in the calculation
        
    Returns:
        str: Path to the written field file
    """
    processor = ScalarFieldProcessor(case_dir)

    # Create field data based on the template
    melt_front_field = template.field.copy()

    print("template field:", melt_front_field.keys())

    # Update header information
    melt_front_field['header']['class'] = 'volScalarField'
    melt_front_field['header']['object'] = template.name
    melt_front_field['header']['location'] = 'postProcessing'
    melt_front_field['header']['dimensions'] = template.dimensions

    # Replace infinite values with a large number for visualization
    viz_data = melt_front_time.copy()
    viz_data[np.isinf(viz_data)] = 1e+30

    # Update the internal field data
    melt_front_field['internal_field'] = viz_data

    # Create output directory
    output_dir = os.path.join(case_dir, 'postProcessing')
    os.makedirs(output_dir, exist_ok=True)
    
    # Write the field file
    output_path = os.path.join('postProcessing', template.name)
    processor.write_field(melt_front_field, output_path)

    print(f"Wrote melt front field to {output_path}")

    for time in times:
        if time == 'template' or time == "0":
            continue

        print(f"\tWriting field for time {time}...")

        output_path = os.path.join(".", str(time), template.name)
        processor.write_field(melt_front_field, output_path)

    return output_path

def calc_stats(melt_front_time):
    """
    Calculate statistics about the fill time.
    
    Args:
        melt_front_time: Array of melt front times
        
    Returns:
        dict: Dictionary of statistics
    """
    # Filter out infinite values (unfilled cells)
    valid_times = melt_front_time[~np.isinf(melt_front_time)]
    
    if len(valid_times) == 0:
        return {
            'min': None,
            'max': None,
            'mean': None,
            'median': None,
            'filled_cells': 0,
            'total_cells': len(melt_front_time),
            'fill_percentage': 0
        }
    
    return {
        'min': np.min(valid_times),
        'max': np.max(valid_times),
        'mean': np.mean(valid_times),
        'median': np.median(valid_times),
        'filled_cells': len(valid_times),
        'total_cells': len(melt_front_time),
        'fill_percentage': len(valid_times) / len(melt_front_time) * 100
    }

def main(case_dir, threshold):
    """Main function to calculate and save melt front times."""
    
    print(f"Processing case: {case_dir}")

    # Extract VoF data from all time steps
    data = extract_data(case_dir)
    
    # Calculate melt front time
    melt_front_time, melt_front_temp = calc_melt_front(data, threshold)

    if melt_front_time is None:
        print("Failed to calculate melt front times. Exiting.")
        return
    
    # Calculate statistics
    stats = calc_stats(melt_front_time)
    
    print("\nFill Time Statistics:")
    print(f"Minimum fill time: {stats['min']:.6f} s")
    print(f"Maximum fill time: {stats['max']:.6f} s")
    print(f"Mean fill time: {stats['mean']:.6f} s")
    print(f"Median fill time: {stats['median']:.6f} s")
    print(f"Filled cells: {stats['filled_cells']} out of {stats['total_cells']} ({stats['fill_percentage']:.2f}%)")
    
    # Write the melt front fields using the template
    write_melt_front_field(case_dir, melt_front_time, data['template'].alpha, data.keys())
    write_melt_front_field(case_dir, melt_front_temp, data['template'].T, data.keys())
    
    print("\nMelt front time calculation completed successfully.")
