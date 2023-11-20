import h5py
import sys
import os
import math
from Utils.I_O import prt, ErrorMsg, ErrorMsgNoQuit
from Utils.I_O import InputHandler

fileName = 'Interface/Modifications.py'

def ModifyTargetFile(targetFile, targetDictPath, targetValue):
  # Retrieve targetFile extension to choose from either Modify* subroutine
  temp      = targetFile.split('.')
  extension = temp[-1]
  # Choose and launch the Modify* subroutine
  if extension == 'h5':
    ModifyH5File(targetFile, targetDictPath, targetValue)
  elif extension == 'input' or extension == 'in':
    ModifyInputFile(targetFile, targetDictPath, targetValue)
  else:
    ErrorMsg('The file ' + targetFile + ' has an unknown extension. Supported extensions are ".in", ".input" and ".h5".', fileName)


def ModifyH5File(targetFile, targetDictPath, targetValue):
  willExit  = False
  # Modifying the sampled mesh file according to the computed sampled value
  try:
    hFile = h5py.File(targetFile, 'r+')
  except:
    ErrorMsgNoQuit('The file ' + targetFile + ' does not exist, or is unreadable as an h5 file.', fileName)
    willExit = True # Necessary because sys.exit is an exception so it is badly handled when inside an "except"
  if willExit:
    sys.exit()
  
  try:
    hFile[targetDictPath][()] = targetValue
    hFile.close()
  except:
    hFile.close()
    ErrorMsgNoQuit('The path ' + targetDictPath + ' does not exit in the file ' + targetFile + '.', fileName)
    willExit = True
  if willExit:
    sys.exit()


def ModifyInputFile(inputFile, dictPath, value):
  # Import input file (agath-like)
  willExit  = False
  try:
    INPUT    = InputHandler(inputFile, modify=True)
  except:
    ErrorMsgNoQuit('The input file ' + inputFile + ' does not exist, or is unreadable as a text file.', fileName)
    willExit = True
  if willExit:
    sys.exit()
  
  # Split the dictPath to get the different keys 
  keyList = dictPath.split('/')
  if len(keyList) < 2:
    ErrorMsgNoQuit('The dictPath ' + dictPath + ' should be at least of length 2 because values are specified in blocks.', fileName)
    willExit = True
  if willExit:
    sys.exit()
  if len(keyList) > 3:
    ErrorMsgNoQuit('The dictPath ' + dictPath + ' should be have a max length of 3 because values are specified in blocks or in blocks/subblocks.', fileName)
    willExit = True
  if willExit:
    sys.exit()
  
  # Navigate the input file and change the value
  INPUT.GetBlock(keyList[0])
  if len(keyList) == 3:
    INPUT.GetSubBlock(keyList[1])
    INPUT.ReplaceFromSpecifier(keyList[2], str(value))
  else:
    INPUT.ReplaceFromSpecifier(keyList[1], str(value))
  INPUT.Write()
  INPUT.Close() # Not necessary ?


def CopyRunFolder(refPath, targetPath):
  if os.path.isdir(refPath):
    os.system("cp -r " + refPath + ' ' + targetPath)
  else:
    ErrorMsg(refPath + ' does not exist or is not a folder.', fileName)


def EraseRunFolder(targetPath):
  os.system("rm -r " + targetPath)
