import os
import hashlib

def __byteStringToHex(byteString):
    return ''.join([hex(ord(c))[2:].zfill(2) for c in byteString])

def randomString(size):
    return __byteStringToHex(os.urandom(size))

def generateActivationToken():
    return __byteStringToHex(os.urandom(25))
def generateSalt():
    return __byteStringToHex(os.urandom(5))

def md5(clear_text):
    return hashlib.sha256(clear_text).hexdigest()

def saltedPassword(clear_password, salt = generateSalt()):
    salted_password = hashlib.sha256(salt + clear_password).hexdigest()
    return salt, salted_password
