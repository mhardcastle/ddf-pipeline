#!usr/bin/python

from astropy.io import fits
import numpy as np
import os
import aplpy
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import pickle
import astropy.units as u
import argparse

PADDING = 0.02
STEP_SIZE = 0.02
BEAM_ATTENUATION_THRESHOLD = 0.3

parser = argparse.ArgumentParser(description="Optimize extract region.")
parser.add_argument("--ra", type=float, default=54.64708333, help="Right Ascension of the target (default: 54.64708333).")
parser.add_argument("--dec", type=float, default=9.965, help="Declination of the target (default: 9.965).")
parser.add_argument("--size", type=float, default=0.4, help="Initial box size in degrees (default: 0.4).")
parser.add_argument("--target_flux", type=float, default=0.2, help="Target total flux in Jy (default: 0.2).")
parser.add_argument("--region_file", type=str, default="region.reg", help="Output DS9 region file (default: region.reg).")
args = parser.parse_args()

ra = args.ra
dec = args.dec
size = args.size
target_flux = args.target_flux
REGION_FILE = args.region_file


def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    """This removes the degenerated dimensions in APLpy 2.X...
    The input must be the object returned by aplpy.FITSFigure().
    `dropaxis` is the index where to start dropping the axis (by default it assumes the 3rd,4th place).
    """
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs


def find_total_quadflux(ra_cen, dec_cen, size_deg, opencat):
    """
    Calculate the total flux of sources within a square region centered at (ra_cen, dec_cen),
    taking into account the curvature of the celestial sphere.

    Parameters:
        ra_cen (float): Center RA of the region in degrees.
        dec_cen (float): Center DEC of the region in degrees.
        size_deg (float): Size of the region in degrees (same for RA and DEC).
        data (astropy Table or structured array): Source catalog containing 'RA', 'DEC', and 'Total_flux'.

    Returns:
        float: Total flux of sources within the region.
    """
    # Convert degrees to radians for calculations
    deg2rad = np.pi / 180.0

    # Adjust RA boundaries based on declination to account for spherical curvature
    ra_offset = size_deg / 2.0 / np.cos(dec_cen * deg2rad)
    ra_min = ra_cen - ra_offset
    ra_max = ra_cen + ra_offset

    # Declination boundaries
    dec_min = dec_cen - size_deg / 2.0
    dec_max = dec_cen + size_deg / 2.0

    # Filter sources within the region
    in_region = (
        (opencat['RA'] >= ra_min) & (opencat['RA'] <= ra_max) &
        (opencat['DEC'] >= dec_min) & (opencat['DEC'] <= dec_max)
    )

    # Calculate the total flux of sources within the region
    total_flux = np.sqrt(np.sum((opencat['Total_flux'][in_region]*opencat['Beam_attenuation'][in_region])**2))

    return total_flux

def find_closest_fields(ra, dec, fieldsdict_path, max_results=10):
    """
    Find the closest entries in the fieldsdict.pkl file to a given RA and DEC.

    Parameters:
        ra (float): Right Ascension of the target in degrees.
        dec (float): Declination of the target in degrees.
        fieldsdict_path (str): Path to the fieldsdict.pkl file.
        max_results (int): Maximum number of closest entries to return.

    Returns:
        list: A list of dictionaries containing 'id', 'ra', 'decl', and 'separation' for the closest entries.
    """
    # Load the fields dictionary from the pickle file
    with open(fieldsdict_path, 'rb') as f:
        fieldsdict = pickle.load(f)

    # Extract RA, DEC, and IDs from the list of dictionaries
    field_ra = np.array([entry['ra'] for entry in fieldsdict])
    field_dec = np.array([entry['decl'] for entry in fieldsdict])
    field_ids = [entry['id'] for entry in fieldsdict]  # Assuming 'id' exists in each entry

    # Create SkyCoord objects for the target and the fields
    target_coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    field_coords = SkyCoord(ra=field_ra * u.degree, dec=field_dec * u.degree, frame='icrs')

    # Calculate separations
    separations = target_coord.separation(field_coords)

    # Sort by separation and get the indices of the closest entries
    closest_indices = np.argsort(separations)[:max_results]

    # Prepare the results
    closest_fields = []
    for idx in closest_indices:
        closest_fields.append({
            'id': field_ids[idx],
            'ra': field_ra[idx],
            'decl': field_dec[idx],
            'separation': separations[idx].degree  # Separation in degrees
        })

    return closest_fields

def get_catalogue_path(field_id):
    """
    Construct the catalogue path for a given field ID.
    """
    return f'../../LoTSS-individual-field-cats-large/Make-matched-catalogues/MASKED_BEAM_CATS/{field_id}_beam_mask.cat.fits'

def download_image(field_id, image_url, output_dir, username, password):
    """
    Download an image using wget with the given credentials.
    """
    image_filename = os.path.join(output_dir, f"{field_id}_image.fits")
    os.system(f"wget {image_url} --user={username} --password={password} -O {image_filename}")
    return image_filename

# Initialize the center of the box
ra_cen, dec_cen = ra, dec  # Start with the initial target RA and DEC
boxsize = size  # Start with the initial box size

# Find the closest fields
closestfields = find_closest_fields(ra, dec, 'fieldsdict.pkl')

# Filter out fields where the average Beam_attenuation is < BEAM_ATTENUATION_THRESHOLD
valid_fields = []
for field in closestfields:
    field_id = field['id']
    catalogue_path = get_catalogue_path(field_id)
    try:
        # Open the catalogue
        hdul = fits.open(catalogue_path)
        catalogue_data = hdul[1].data  # Assume the first extension contains the data

        # Calculate the average Beam_attenuation for the initial box in this field
        deg2rad = np.pi / 180.0
        ra_offset = size / 2.0 / np.cos(dec * deg2rad)
        ra_min = ra - ra_offset
        ra_max = ra + ra_offset
        dec_min = dec - size / 2.0
        dec_max = dec + size / 2.0

        # Filter sources within the region
        in_region = (
            (catalogue_data['RA'] >= ra_min) & (catalogue_data['RA'] <= ra_max) &
            (catalogue_data['DEC'] >= dec_min) & (catalogue_data['DEC'] <= dec_max)
        )

        # Calculate the average Beam_attenuation for sources in the region
        if np.any(in_region):  # Ensure there are sources in the region
            avg_beam_attenuation = np.mean(catalogue_data['Beam_attenuation'][in_region])
        else:
            avg_beam_attenuation = 0.0  # No sources in the region

        # Check if the average Beam_attenuation meets the threshold
        if avg_beam_attenuation >= BEAM_ATTENUATION_THRESHOLD:
            valid_fields.append(field)  # Keep the field if the average Beam_attenuation is valid
        else:
            print(f"Field {field_id} removed due to low average Beam_attenuation: {avg_beam_attenuation:.4f}")

    except FileNotFoundError:
        print(f"Catalogue not found for field ID: {field_id}")
    except Exception as e:
        print(f"Error processing field ID {field_id}: {e}")

# Update closestfields to only include valid fields
closestfields = valid_fields

# Filter out fields that return NaN flux
valid_fields = []
for field in closestfields:
    field_id = field['id']
    catalogue_path = get_catalogue_path(field_id)
    try:
        # Open the catalogue
        hdul = fits.open(catalogue_path)
        catalogue_data = hdul[1].data  # Assume the first extension contains the data

        # Calculate the flux for the initial box in this field
        field_flux = find_total_quadflux(ra, dec, size, catalogue_data)

        # Check if the flux is NaN
        if not np.isnan(field_flux):
            valid_fields.append(field)  # Keep the field if the flux is valid
        else:
            print(f"Field {field_id} removed due to NaN flux.")

    except FileNotFoundError:
        print(f"Catalogue not found for field ID: {field_id}")
    except Exception as e:
        print(f"Error processing field ID {field_id}: {e}")

# Update closestfields to only include valid fields
closestfields = valid_fields

# Load all catalogues into memory before the loop
catalogues = {}
for field in closestfields:
    field_id = field['id']
    catalogue_path = get_catalogue_path(field_id)   
    try:
        # Open the catalogue and store it in memory
        hdul = fits.open(catalogue_path)
        catalogue_data = hdul[1].data  # Assume the first extension contains the data
        catalogues[field_id] = catalogue_data
        hdul.close()  # Close the file after loading the data
    except FileNotFoundError:
        print(f"Catalogue not found for field ID: {field_id}")
    except Exception as e:
        print(f"Error loading catalogue for field ID {field_id}: {e}")

# Ensure the candidate box always contains the original region
while True:
    all_fields_satisfied = True  # Flag to check if all fields meet the target flux
    best_ra_cen, best_dec_cen = ra_cen, dec_cen  # Track the best center for this iteration
    max_flux_per_field = {}  # Dictionary to store the maximum flux for each field

    # Define the RA and DEC ranges for the current box
    rarange = np.arange(ra_cen - (boxsize / 2.0) + PADDING, ra_cen + (boxsize / 2.0) - PADDING, STEP_SIZE)
    decrange = np.arange(dec_cen - (boxsize / 2.0) + PADDING, dec_cen + (boxsize / 2.0) - PADDING, STEP_SIZE)

    # Initialize the progress bar
    total_iterations = len(rarange) * len(decrange)
    progress_bar = tqdm(total=total_iterations, desc=f"Box Size: {boxsize:.4f} degrees", unit="iteration")

    # Loop over candidate centers within the current box
    for ra_candidate in rarange:
        for dec_candidate in decrange:
            # Update the progress bar
            progress_bar.update(1)

            # Adjust the candidate box to ensure it fully contains the original region
            ra_min_candidate = min(ra_candidate - boxsize / 2.0, ra - size / 2.0)
            ra_max_candidate = max(ra_candidate + boxsize / 2.0, ra + size / 2.0)
            dec_min_candidate = min(dec_candidate - boxsize / 2.0, dec - size / 2.0)
            dec_max_candidate = max(dec_candidate + boxsize / 2.0, dec + size / 2.0)

            # Recalculate the candidate box center and size
            adjusted_ra_cen = (ra_min_candidate + ra_max_candidate) / 2.0
            adjusted_dec_cen = (dec_min_candidate + dec_max_candidate) / 2.0
            adjusted_boxsize = max(ra_max_candidate - ra_min_candidate, dec_max_candidate - dec_min_candidate)

            # Check if the flux in each field satisfies the target flux
            all_fields_satisfied = True  # Reset for this candidate center
            fields_to_remove = []  # Track fields to remove due to NaN flux
            for field in closestfields:
                field_id = field['id']

                try:
                    # Retrieve the catalogue data from memory
                    catalogue_data = catalogues[field_id]

                    # Calculate the flux for the adjusted box in this field
                    field_flux = find_total_quadflux(adjusted_ra_cen, adjusted_dec_cen, adjusted_boxsize, catalogue_data)

                    # Check if the flux is NaN
                    if np.isnan(field_flux):
                        print(f"Field ID: {field_id} removed due to NaN flux.")
                        fields_to_remove.append(field)
                        all_fields_satisfied = False
                        break  # Exit the loop for this candidate center

                    # Update the maximum flux for this field
                    if field_id not in max_flux_per_field or field_flux > max_flux_per_field[field_id]:
                        max_flux_per_field[field_id] = field_flux

                    # Check if the flux in this field meets the target
                    if field_flux < target_flux:
                        all_fields_satisfied = False  # This field does not meet the target
                        break  # No need to check further fields for this candidate center

                except KeyError:
                    print(f"Catalogue data not found in memory for field ID: {field_id}")
                    fields_to_remove.append(field)
                    all_fields_satisfied = False
                    break
                except Exception as e:
                    print(f"Error processing field ID {field_id}: {e}")
                    fields_to_remove.append(field)
                    all_fields_satisfied = False
                    break

            # Remove fields with NaN flux from closestfields
            for field in fields_to_remove:
                closestfields.remove(field)

            # If all fields satisfy the target flux, update the best center
            if all_fields_satisfied:
                print(f"Final box parameters before exiting loop:")
                print(f"  Center RA: {adjusted_ra_cen:.6f}, DEC: {adjusted_dec_cen:.6f}")
                print(f"  Box Size: {adjusted_boxsize:.4f} degrees")
                break  # Exit the loop since we found a valid center

        if all_fields_satisfied:
            break  # Exit the outer loop as well

    # Close the progress bar
    progress_bar.close()

    # If all fields satisfy the target flux, stop growing the box
    if all_fields_satisfied:
        break

    # Print the maximum flux found in each field for the current box size
    print(f"Growing box: New Box Size = {boxsize:.4f} degrees")
    for field_id, max_flux in max_flux_per_field.items():
        print(f"  Field ID: {field_id}, Maximum Flux: {max_flux:.4f} Jy")

    # Grow the box size if the target flux is not yet reached in all fields
    boxsize *= 1.1  # Increase the box size by 10%

# Ensure the final box contains the original box
final_boxsize = max(adjusted_boxsize * 1.0, size)  # Ensure the final box size is at least as large as the original box
final_ra_cen = adjusted_ra_cen
final_dec_cen = adjusted_dec_cen

# Refine the final box center to maximize the flux while ensuring it contains the original box
print("\nRefining the final box center to maximize flux while ensuring it contains the original box...")
best_flux = -np.inf  # Initialize with a very low value
best_ra_cen = final_ra_cen
best_dec_cen = final_dec_cen

# Define the RA and DEC ranges for the final box
rarange = np.arange(final_ra_cen - (final_boxsize / 2.0) + PADDING, final_ra_cen + (final_boxsize / 2.0) - PADDING, STEP_SIZE)
decrange = np.arange(final_dec_cen - (final_boxsize / 2.0) + PADDING, final_dec_cen + (final_boxsize / 2.0) - PADDING, STEP_SIZE)

# Iterate through the RA and DEC ranges to find the center with the maximum flux
for ra_candidate in rarange:
    for dec_candidate in decrange:
        # Ensure the candidate box contains the original box
        ra_min_candidate = ra_candidate - final_boxsize / 2.0
        ra_max_candidate = ra_candidate + final_boxsize / 2.0
        dec_min_candidate = dec_candidate - final_boxsize / 2.0
        dec_max_candidate = dec_candidate + final_boxsize / 2.0

        if (
            ra_min_candidate <= ra - size / 2.0 and
            ra_max_candidate >= ra + size / 2.0 and
            dec_min_candidate <= dec - size / 2.0 and
            dec_max_candidate >= dec + size / 2.0
        ):
            total_flux = 0.0  # Initialize total flux for this candidate center
            for field in closestfields:
                field_id = field['id']
                try:
                    # Retrieve the catalogue data from memory
                    catalogue_data = catalogues[field_id]

                    # Calculate the flux for the candidate box in this field
                    field_flux = find_total_quadflux(ra_candidate, dec_candidate, final_boxsize, catalogue_data)

                    # Add the field flux to the total flux
                    total_flux += field_flux

                except KeyError:
                    print(f"Catalogue data not found in memory for field ID: {field_id}")
                except Exception as e:
                    print(f"Error processing field ID {field_id}: {e}")

            # Update the best center if this candidate has a higher total flux
            if total_flux > best_flux:
                best_flux = total_flux
                best_ra_cen = ra_candidate
                best_dec_cen = dec_candidate

# Update the final center with the best center found
final_ra_cen = best_ra_cen
final_dec_cen = best_dec_cen
final_boxsize = adjusted_boxsize *1.05 # Stretch the box a tiny bit to try to prevent sources being cut in half.

# Print final results
print(f"\nFinal box size (adjusted): {final_boxsize:.4f} degrees")
print(f"Final center (refined): RA = {final_ra_cen:.6f}, DEC = {final_dec_cen:.6f}")
print(f"Maximum flux contained: {best_flux:.4f} Jy")

# Calculate and print the flux in each field for the final box
print("\nFlux in each field for the final box:")
for field in closestfields:
    field_id = field['id']

    try:
        # Retrieve the catalogue data from memory
        catalogue_data = catalogues[field_id]

        # Calculate the flux for the final box in this field
        field_flux = find_total_quadflux(final_ra_cen, final_dec_cen, final_boxsize, catalogue_data)

        # Print the flux for this field
        print(f"Field ID: {field_id}, Flux: {field_flux:.4f} Jy")

    except KeyError:
        print(f"Catalogue data not found in memory for field ID: {field_id}")
    except Exception as e:
        print(f"Error processing field ID {field_id}: {e}")

# Write the final region to a DS9 region file
with open(REGION_FILE, 'w') as g:
    g.write('# Region file format: DS9 version 4.1\n')
    g.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    g.write('fk5\n')
    # Write the final box region with the adjusted size
    g.write(f'box({final_ra_cen},{final_dec_cen},{final_boxsize*3600:.2f}",{final_boxsize*3600:.2f}",0.0) # color=magenta\n')
    # Write the original target point
    g.write(f'point({ra},{dec}) # point = circle color=red\n')
    # Write the initial box region with the original size
    g.write(f'box({ra},{dec},{size*3600:.2f}",{size*3600:.2f}",0.0) # color=green\n')


print(final_boxsize, final_ra_cen, final_dec_cen)

# Read the password from ~/.surveys file
surveys_credentials_path = os.path.expanduser("~/.surveys")
if not os.path.exists(surveys_credentials_path):
    raise FileNotFoundError("The ~/.surveys file does not exist. Please create it with your password.")

with open(surveys_credentials_path, 'r') as f:
    surveys_password = f.readline().strip()

# The username is always 'surveys'
surveys_username = "surveys"

# Directory to store the downloaded images
output_dir = "downloaded_images"
os.makedirs(output_dir, exist_ok=True)

# Loop through the remaining fields in closestfields
for field in closestfields:
    field_id = field['id']
    image_url = f"https://lofar-surveys.org/downloads/DR3/fields/{field_id}/image_full_ampphase_di_m.NS.app.restored.fits"
    png_filename = os.path.join(output_dir, f"{field_id}_postage_stamp.png")

    # Download the image using the new function
    print(f"Downloading image for field {field_id}...")
    image_filename = download_image(field_id, image_url, output_dir, surveys_username, surveys_password)

    # Check if the file was downloaded successfully
    if not os.path.exists(image_filename):
        print(f"Failed to download image for field {field_id}. Skipping...")
        continue

    # Create a postage stamp image using APLpy
    print(f"Creating postage stamp for field {field_id}...")
    try:
        fig = aplpy.FITSFigure(image_filename)
        fig.show_grayscale(stretch='linear', invert=True)

        # Use the final region file directly
        fig.show_regions(REGION_FILE)

        # Save the postage stamp image as a PNG
        fig.save(png_filename)
        fig.close()
        print(f"Postage stamp saved: {png_filename}")

    except Exception as e:
        print(f"Error creating postage stamp for field {field_id}: {e}")

