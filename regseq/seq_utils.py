def stitch(s):
    """Combine the mutated sequence with barcode."""
    return s['seq'] + s['tag']


def format_string(x):
    'basic function to format output string'
    return '%10.6f' %x
