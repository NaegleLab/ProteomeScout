import codecs
import chardet
import shutil

def recode(infilename, outfilename, source_encoding, target_encoding):
    BLOCKSIZE = 1048576 # or some other, desired size in bytes
    with codecs.open(infilename, 'r', source_encoding) as sourceFile:
        with codecs.open(outfilename, 'w', target_encoding) as targetFile:
            while True:
                contents = sourceFile.read(BLOCKSIZE)
                if not contents:
                    break
                targetFile.write(contents)

def test_file_encoding(filename, encoding):
    try:
        with codecs.open(filename, 'r', encoding) as f:
            for _line in f: pass
    except UnicodeDecodeError:
        return False
    return True

def convert_encoding_to_utf8(infilename, outfilename):
    prediction = chardet.detect(infilename)

    valid_encodings = set(['ascii', 'utf-8', 'utf8'])
    predicted_encoding = prediction['encoding'].lower()

    if predicted_encoding not in valid_encodings:
        recode(infilename, outfilename, predicted_encoding, 'utf8')
    else:
        i = 0
        alternate_encodings = ['utf8','cp1252','utf16','utf32']
        while not test_file_encoding(infilename, alternate_encodings[i]):
            i+=1

        if alternate_encodings[i] == 'utf8':
            shutil.copyfile(infilename, outfilename)
        else:
            recode(infilename, outfilename, alternate_encodings[i], 'utf8')
