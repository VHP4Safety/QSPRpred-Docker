import zlib

def encode(string: str, encoding='utf-8'):
    """
    Compresses a string
    Arguments:
        string (str)              : string to compress  
    Returns:
        compressed (str)          : compressed string
    """ 
    if encoding == 'utf-8':
        return zlib.compress(string.encode('utf-8')).hex()

def decode(string: str, encoding='utf-8'):
    """
    Decompresses a compressed string
    Arguments:
        compressed (str)          : compressed string
    Returns:
        string (str)              : decompressed string
    """
    if encoding == 'utf-8':
        return zlib.decompress(bytes.fromhex(string)).decode('utf-8')
