import re

def split_filenames(filenames):
    filenames = filenames.lstrip("[").rstrip("]").split(", ")
    return filenames

def map_chromosome_filename(filename):
    return re.search(r"(?:chr)(\d{1,2}|X|Y)", filename).groups()[0]