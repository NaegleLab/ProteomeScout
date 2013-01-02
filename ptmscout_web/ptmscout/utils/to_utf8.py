import codecs
import chardet
import shutil


def convert_encoding_to_utf8(infilename, outfilename):
    prediction = chardet.detect(infilename)

    valid_encodings = set(['ascii', 'utf-8'])
    predicted_encoding = prediction['encoding'].lower()

    BLOCKSIZE = 1048576 # or some other, desired size in bytes

    if predicted_encoding not in valid_encodings:
        with codecs.open(infilename, 'r', predicted_encoding) as sourceFile:
            with codecs.open(outfilename, 'w', 'utf-8') as targetFile:
                while True:
                    contents = sourceFile.read(BLOCKSIZE)
                    if not contents:
                        break
                    targetFile.write(contents)
    else:
        shutil.copyfile(infilename, outfilename)
