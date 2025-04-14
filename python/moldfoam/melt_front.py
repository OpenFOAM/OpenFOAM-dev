import os
import numpy as np
from pathlib import Path
from moldfoam.scalar import ScalarFieldProcessor

def extract_alpha_data(case_dir):
    """
    Extract alpha field data from all time directories in the case.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        
    Returns:
        dict: Dictionary mapping time values to alpha field data
    """
    case_path = Path(case_dir)
    processor = ScalarFieldProcessor(case_dir)
    
    # Get time directories sorted numerically
    time_dirs = sorted([d for d in os.listdir(case_dir) 
                     if d.replace('.', '', 1).isdigit()], 
                     key=lambda x: float(x))
    
    # Dictionary to store alpha values by time
    alpha_by_time = {}
    
    # Read alpha field for each time step
    print(f"Reading alpha field data from {len(time_dirs)} time directories...")
    for time_dir in time_dirs:
        time_val = float(time_dir)
    
        alpha_file = Path(time_dir) / "alpha.melt"
        try:
            # Read the alpha field
            alpha_field = processor.read_field(str(alpha_file))
            
            # Store only the internal field data for efficiency
            alpha_by_time[time_val] = {
                'data': alpha_field['internal_field'],
                'path': str(alpha_file)
            }
            
            if len(alpha_by_time) == 1:
                # Store the first alpha field completely to use as a template later
                alpha_by_time['template'] = alpha_field
        except Exception as e:
            print(f"Error reading {alpha_file}: {str(e)}")
    
    print(f"Successfully read alpha data from {len(alpha_by_time)-1} time steps")
    return alpha_by_time


def cross_time(alpha_curr, alpha_prev, t_curr, t_prev):
    alpha_slope = (alpha_curr - alpha_prev) / (t_curr - t_prev)
    delta_t = (0.5 - alpha_prev) / alpha_slope
    return t_prev + delta_t


def calc_melt_front_time(alpha_by_time):
    """
    Calculate the time when alpha crosses 0.5 for each cell.
    
    Args:
        alpha_by_time: Dictionary mapping time values to alpha field data
        
    Returns:
        numpy.ndarray: Array of melt front times for each cell
    """
    # Remove the template entry before sorting times
    alpha_data = {k: v for k, v in alpha_by_time.items() if k != 'template'}
    
    # Get sorted time steps
    times = sorted(alpha_data.keys())
    
    if not times:
        raise ValueError("No time steps found with valid alpha data")
    
    # Get number of cells from the first time step
    first_alpha = alpha_data[times[1]]['data']
    if isinstance(first_alpha, (int, float)):
        print("Warning: Alpha field is uniform. Cannot calculate melt front times.")
        return None
    
    num_cells = len(first_alpha)
    print(f"Calculating melt front times for {num_cells} cells...")
    
    # Initialize melt_front_time with a large value
    melt_front_time = np.ones(num_cells) * float('inf')
    
    # For each time step
    for i in range(1, len(times)):
        t_curr = times[i]
        t_prev = times[i-1]
        alpha_curr = alpha_data[t_curr]['data']
    
        if i == 1:
            alpha_prev = np.zeros_like(alpha_curr)
        else:
            alpha_prev = alpha_data[t_prev]['data']
        
        # Find cells where alpha crosses 0.5 between these time steps
        # We need the cells where alpha was below 0.5 and now is at or above 0.5
        crossings = np.where((alpha_prev < 0.5) & (alpha_curr >= 0.5))[0]
    
        for cell in crossings:
            # Linear interpolation to find more precise crossing time
            crossing_time = cross_time(
                alpha_curr[cell], alpha_prev[cell], t_curr, t_prev
            )
            
            # Update melt_front_time if this is the first crossing
            if melt_front_time[cell] == float('inf'):
                melt_front_time[cell] = crossing_time
    
    # Count cells that were filled
    filled_cells = np.sum(melt_front_time < float('inf'))
    print(f"Found melt front times for {filled_cells} out of {num_cells} cells")

    melt_front_time[np.isinf(melt_front_time)] = times[-1]
    
    return melt_front_time

def write_melt_front_field(case_dir, melt_front_time, alpha_template):
    """
    Write the melt front time field to an OpenFOAM field file.
    
    Args:
        case_dir: Path to the OpenFOAM case directory
        melt_front_time: Array of melt front times
        alpha_template: Template alpha field to use for metadata
        
    Returns:
        str: Path to the written field file
    """
    processor = ScalarFieldProcessor(case_dir)
    
    # Create field data based on the template
    melt_front_field = alpha_template.copy()
    
    # Update header information
    melt_front_field['header']['class'] = 'volScalarField'
    melt_front_field['header']['object'] = 'meltFrontTime'
    melt_front_field['header']['location'] = 'postProcessing'
    # Set dimensions to [0 0 1 0 0 0 0] for time in seconds
    melt_front_field['header']['dimensions'] = [0, 0, 1, 0, 0, 0, 0]
    
    # Replace infinite values with a large number for visualization
    viz_data = melt_front_time.copy()
    viz_data[np.isinf(viz_data)] = 1e+30
    
    # Update the internal field data
    melt_front_field['internal_field'] = viz_data
    
    # Create output directory
    output_dir = os.path.join(case_dir, 'postProcessing')
    os.makedirs(output_dir, exist_ok=True)
    
    # Write the field file
    output_path = os.path.join('postProcessing', 'meltFrontTime')
    processor.write_field(melt_front_field, output_path)
    
    print(f"Wrote melt front time field to {output_path}")

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

def main(case_dir):
    """Main function to calculate and save melt front times."""
    
    print(f"Processing case: {case_dir}")
    
    # Extract VoF data from all time steps
    data = extract_alpha_data(case_dir)
    
    # Calculate melt front time
    melt_front_time = calc_melt_front_time(data)
    
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
    
    # Write the melt front time field using the template
    write_melt_front_field(case_dir, melt_front_time, data['template'])
    
    print("\nMelt front time calculation completed successfully.")
