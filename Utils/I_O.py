import numpy as np
import sys
import copy
import inspect
from termcolor import colored as cld


###################################################
# Colorful terminal outputs with verbosity option #
###################################################
def prt(message, colour, v):
    if v >= 0:
        if colour == 'None':
            print(message)
        else:
            print(cld(message, colour))


def ErrorMsg(message, fileName):
    routineName = inspect.currentframe().f_back.f_code.co_name
    prt('', 'red', True)
    prt('Location: ' + fileName + '/' + routineName, 'red', True)
    prt('Error: ' + message, 'red', True)
    prt('Exiting...', 'red', True)
    sys.exit()


def ErrorMsgNoQuit(message, fileName):
    routineName = inspect.currentframe().f_back.f_code.co_name
    prt('', 'red', True)
    prt('Location: ' + fileName + '/' + routineName, 'red', True)
    prt('Error: ' + message, 'red', True)
    prt('Exiting...', 'red', True)


def WarningMsg(message, fileName, verbosity):
    routineName = inspect.currentframe().f_back.f_code.co_name
    prt('', 'yellow', verbosity)
    prt('Location: ' + fileName + '/' + routineName, 'yellow', verbosity)
    prt('Warning: ' + message, 'yellow', verbosity)


############################################
#              Load CSV file               #
############################################
def LoadCSV(fileName, delimiter=',', nHeaderLines=6):
    X = []
    Y = []

    with open(fileName, 'r') as f:
        content = f.readlines()
        for lineIdx in range(len(content)):
            if lineIdx >= nHeaderLines:
                line = content[lineIdx]
                if delimiter == " ":
                    row = line.split()
                else:
                    row = line.split(delimiter)
                X.append(float(row[0]))
                Y.append(float(row[1]))
    return np.array(X), np.array(Y)


####################################################################
# Input file handler module from Matteo Gelain (EM2C/ONERA/SAFRAN) #
# Suitable for Agath-like files                                    #
# Shared on October, 21st 2019                                     #
####################################################################
class AgathInputHandler:
    # TODO: Make the "keyword" file optional and if it is not present, guess the type of the data
    # TODO: Increase genericity by proposing arguments to change the comment sign (!, #, ...) and the separator between the key and the value(s) (:, =, ...) => Choose from Agath or AVBP input file types at initialization
    def __init__(self, filename, **kwargs):
        self.filename = filename
        if 'modify' in kwargs and kwargs['modify'] == True:
            self.file = open(filename, 'r+')
            self.Modify = True
        else:
            self.file = open(filename, 'r')
            self.Modify = False
        self.all_lines_raw = self.file.readlines()
        self.all_lines = copy.deepcopy(self.all_lines_raw)
        self.all_lines = self.CleanLines(self.all_lines)
        self.currentLines = self.all_lines

    def Close(self):
        self.file.close()
        for key in self.__dict__.keys():
            exec('del self.%s' % key)

    def GetBlock(self, blockName):
        self.blockName = blockName
        self.blockLines = []
        start_index = - 1
        end_index = - 1
        for i in range(len(self.all_lines)):
            line = self.all_lines[i]
            if 'BLOCK' in line and 'END' not in line and blockName in line and 'SUBBLOCK' not in line:
                start_index = i
            elif 'BLOCK' in line and 'END' in line and 'SUBBLOCK' not in line:
                end_index = i
                if start_index >= 0:
                    break

        if start_index < 0 and end_index < 0:
            prt('BLOCK ' + blockName + ' not found.', 'red', True)
            prt('Exiting.', 'red', 1)
            sys.exit()
        elif start_index < 0:
            print('Beginning of BLOCK ' + blockName + ' not found.')
            print('Stopping.')
            sys.exit()
        elif end_index < 0:
            print('End of BLOCK ' + blockName + ' not found.')
            print('Stopping.')
            sys.exit()
        else:
            self.blockLines = self.all_lines[start_index + 1:end_index]
            self.currentLines = self.blockLines
            self.currentBlock = self.blockName
            if 'subBlockName' in self.__dict__.keys():
                del self.subBlockName
                del self.subBlockLines
        return

    def GetSubBlock(self, blockName):
        self.subBlockName = blockName
        self.subBlockLines = []
        if 'blockLines' not in self.__dict__.keys():
            print('A BLOCK of inputs is needed before looking for a SUBBLOCK.')
            print('Stopping.')
            sys.exit()
        elif self.blockLines == [] and 'blockName' in self.__dict__.keys():
            print('It seems that no lines were found for BLOCK: ' + self.blockName)
            print('I cannot look for SUBBLOCK ' + blockName + '.')
            print('Stopping.')
            sys.exit()
        start_index = - 1
        end_index = - 1
        for i in range(len(self.blockLines)):
            line = self.blockLines[i]
            if 'SUBBLOCK' in line and 'END' not in line and blockName in line:
                start_index = i
            elif 'SUBBLOCK' in line and 'END' in line:
                end_index = i
                if start_index >= 0:
                    break

        if start_index < 0 and end_index < 0:
            print('SUBBLOCK ' + blockName + ' not found.')
            print('Stopping.')
            sys.exit()
        elif start_index < 0:
            print('Beginning of SUBBLOCK ' + blockName + ' not found.')
            print('Stopping.')
            sys.exit()
        elif end_index < 0:
            print('End of SUBBLOCK ' + blockName + ' not found.')
            print('Stopping.')
            sys.exit()
        else:
            self.subBlockLines = self.blockLines[start_index + 1:end_index]
            self.currentLines = self.subBlockLines
            self.currentBlock = self.subBlockName

    def GetValueFromKey(self, keyName, typ):
        self.keyName = keyName
        self.Value = []
        if 'currentBlock' not in self.__dict__.keys():
            print('A block of inputs is needed before looking for a value.')
            print('Stopping.')
            sys.exit()
        elif self.currentLines == [] and 'currentBlock' in self.__dict__.keys():
            print('It seems that no lines were found for block : ' + self.currentBlock)
            print('I cannot look for ' + keyName)
            print('Stopping.')
            sys.exit()
        else:
            found = False
            for line in self.currentLines:
                if 'SUBBLOCK' not in line:
                    line = line.strip()
                    line = line.split('=')
                    if len(line) < 2:
                        print('Problem for key ' + keyName + ' in block ' + self.currentBlock + '.')
                        print('A little bird told me that you might have to use "=" as separator...')
                        sys.exit()
                    if line[0].strip() == keyName:
                        found = True
                        del line[0]
                        line[-1] = line[-1].strip()
                        line = line[-1].split(' ')
                        if line[0] != '-->':
                            for word in line:
                                word.strip()
                                if typ == 'char':
                                    Value = word
                                    if word == '[]':
                                        Value = []
                                elif typ == 'int':
                                    try:
                                        Value = int(word)
                                    except:
                                        if word == '[]':
                                            Value = []
                                        else:
                                            print(
                                                'Input for key ' + keyName + ' in block ' + self.currentBlock + ' is not an integer as specified.')
                                            print('For an empty list, please write "[]".')
                                            print('Stopping.')
                                            sys.exit()
                                elif typ == 'real':
                                    try:
                                        Value = float(word)
                                    except:
                                        if word == '[]':
                                            Value = []
                                        else:
                                            print(
                                                'Input for key ' + keyName + ' in block ' + self.currentBlock + ' is not a real as specified.')
                                            print('For an empty list, please write "[]".')
                                            print('Stopping.')
                                            sys.exit()
                                elif typ == 'logical':
                                    if word == 'True':
                                        Value = True
                                    elif word == 'False':
                                        Value = False
                                    else:
                                        if word == '[]':
                                            Value = []
                                        else:
                                            print(
                                                'Input for key ' + keyName + ' in block ' + self.currentBlock + ' should be either True or False.')
                                            print('For an empty list, please write "[]".')
                                            print('Got ' + word + ' instead.')
                                            print('Stopping.')
                                            sys.exit()
                                else:
                                    print(
                                        'Type for key ' + keyName + ' in block ' + self.currentBlock + ' seems not to be recognised.')
                                    print('Type should be in [char, int, real, logical]. You entered : ' + str(typ))
                                    print('Stopping.')
                                    sys.exit()

                                if len(line) > 1:
                                    self.Value.append(Value)
                                elif len(line) == 1:
                                    self.Value = Value
                                else:
                                    print(
                                        'No values for key ' + keyName + ' in block ' + self.currentBlock + ' were found.')
                                    print('Stopping.')
                                    sys.exit()
                            return self.Value

                        else:  # Handles references to other specifiers
                            # If several references are given with : "* = --> ref1 ref2", a list of the values is formed
                            # The path must be complete, i.e. be of the form blockID/(subBlockID/)specifier
                            if len(line) < 2:
                                prt('', 'red', True)
                                prt('Utils/I_O.py/InputHandler/GetValueFromKey', 'red', True)
                                prt('Error: No ref given', 'red', True)
                                prt('Exiting', 'red', True)
                                sys.exit()
                            else:
                                blockName_old = self.blockName
                                subBlockName_old = ''
                                if 'subBlockName' in self.__dict__.keys():
                                    subBlockName_old = self.subBlockName

                                if len(line) == 2:  # If there is only one reference
                                    path = line[1].split('/')
                                    if len(path) == 2:
                                        blockName = path[0]
                                        keyName = path[1]
                                        self.GetBlock(blockName)
                                        self.GetValueFromKey(keyName, typ)
                                    elif len(path) == 3:
                                        blockName = path[0]
                                        subBlockName = path[1]
                                        keyName = path[2]
                                        self.GetBlock(blockName)
                                        self.GetValueFromKey(keyName, typ)
                                    else:
                                        prt('', 'red', True)
                                        prt('Utils/I_O.py/InputHandler/GetValueFromKey', 'red', True)
                                        prt('Error: ' + line[
                                            1] + ' is not a valid reference. Valid format is "blockID/(subBlockID/)specifier"',
                                            'red', True)
                                        prt('Exiting', 'red', True)
                                        sys.exit()
                                else:
                                    myValue = []
                                    for refIdx in range(len(line) - 1):
                                        path = line[refIdx + 1].split('/')
                                        if len(path) == 2:
                                            blockName = path[0]
                                            keyName = path[1]
                                            self.GetBlock(blockName)
                                            self.GetValueFromKey(keyName, typ)
                                            myValue.append(self.Value)
                                        elif len(path) == 3:
                                            blockName = path[0]
                                            subBlockName = path[1]
                                            keyName = path[2]
                                            self.GetBlock(blockName)
                                            self.GetSubBlock(subBlockName)
                                            self.GetValueFromKey(keyName, typ)
                                            myValue.append(self.Value)
                                        else:
                                            prt('', 'red', True)
                                            prt('Utils/I_O.py/InputHandler/GetValueFromKey', 'red', True)
                                            prt('Error: ' + line[
                                                1] + ' is not a valid reference. Valid format is "blockID/(subBlockID/)specifier"',
                                                'red', True)
                                            prt('Exiting', 'red', True)
                                            sys.exit()
                                    self.Value = myValue

                            self.GetBlock(blockName_old)
                            if subBlockName_old != '':
                                self.GetSubBlock(subBlockName_old)
                            return self.Value

            if found == False:
                print('Key ' + keyName + ' not found in BLOCK ' + self.currentBlock)
                print('Stopping.')
                sys.exit()

    def GetVariableFromKey(self, KEYWORDS, keyName):
        typ = KEYWORDS.GetCharacterFromKey(keyName)
        if isinstance(typ, list):
            typ = typ[0]
        Value = self.GetValueFromKey(keyName, typ)
        return Value

    def GetCharacterFromKey(self, keyName):
        Value = self.GetValueFromKey(keyName, 'char')
        return Value

    def GetIntegerFromKey(self, keyName):
        Value = self.GetValueFromKey(keyName, 'int')
        return Value

    def GetRealFromKey(self, keyName):
        Value = self.GetValueFromKey(keyName, 'real')
        return Value

    def GetLogicalFromKey(self, keyName):
        Value = self.GetValueFromKey(keyName, 'logical')
        return Value

    def GetKeysForBlock(self, **kwargs):
        blockKeys = []
        if 'currentBlock' not in self.__dict__.keys():
            print('A block of inputs is needed before looking for block keys.')
            print('Stopping.')
            sys.exit()
        elif self.currentLines == [] and 'currentBlock' in self.__dict__.keys():
            print('It seems that no lines were found for block : ' + self.currentBlock)
            print('I cannot look for keys.')
            print('Stopping.')
            sys.exit()
        else:
            for line in self.currentLines:
                line = line.strip()
                if line != '' and 'SUBBLOCK' not in line:
                    if 'key' in kwargs.keys():
                        if kwargs['key'] in line:
                            line = line.split('=')
                            if len(line) < 2:
                                print('Problem for block ' + self.currentBlock + ' while looking for keys.')
                                print('A little bird told me that you might have to use "=" as separator...')
                                sys.exit()
                            blockKeys.append(line[0].strip())
                    else:
                        line = line.split('=')
                        if len(line) < 2:
                            print('Problem for block ' + self.currentBlock + ' while looking for keys.')
                            print('A little bird told me that you might have to use "=" as separator...')
                            sys.exit()
                        blockKeys.append(line[0].strip())
        self.blockKeys = blockKeys
        return blockKeys

    def CleanLines(self, lines):
        # TODO: Choose AVBP or Agath version when initializing the handler
        index = []
        for i in range(len(lines)):
            lines[i] = lines[i].strip()
            if lines[i] == '' or lines[i].startswith('#'):  # Agath version
                # if lines[i] == '' or lines[i].startswith('!'): # AVBP version
                index.append(i)
        lines_temp = [lines[i] for i in range(len(lines)) if i not in index]
        lines = lines_temp
        for i in range(len(lines)):
            lines[i] = lines[i].split('#')[0]  # Agath version
        # lines[i] = lines[i].split('!')[0] # AVBP version
        return lines

    def GetSubBlockNames(self, blockName=''):
        subBlockNames = []
        self.subBlockName = blockName
        self.subBlockLines = []
        if 'blockLines' not in self.__dict__.keys():
            print('A block of inputs is needed before looking for a SUBBLOCK names.')
            print('Stopping.')
            sys.exit()
        elif self.blockLines == [] and 'blockName' in self.__dict__.keys():
            print('It seems that no lines were found for block : ' + self.blockName)
            print('I cannot look for SUBBLOCKs.')
            print('Stopping.')
            sys.exit()
        else:
            for i in range(len(self.blockLines)):
                line = self.blockLines[i]
                line = line.strip()
                if 'SUBBLOCK' in line and 'END' not in line and blockName in line:
                    line = line.split('SUBBLOCK')
                    subBlockNames.append(line[1].strip())
        self.subBlockNames = subBlockNames
        return subBlockNames

    def ReplaceFromKey(self, key, value):
        value = str(value)
        if self.Modify == True:
            for i in range(len(self.all_lines_raw)):
                line = self.all_lines_raw[i]
                if key in line:
                    self.all_lines_raw[i] = self.all_lines_raw[i].replace(key, value)
            # sys.stdout.write(line)
        else:
            print('Function ReplaceFromKey can only be used in modify mode.')
            print('Please add "modify=True" while opening file.')
            print('Stopping.')
            sys.exit()
        return

    def ReplaceFromSpecifier(self, key, value, **kwargs):
        check = False
        check_block = True
        check_subBlock = True
        if 'blockName' in self.__dict__.keys():
            check_block = False
            if 'subBlockName' in self.__dict__.keys():
                check_subBlock = False

        if isinstance(value, list):
            string = ''
            for val in value:
                string = string + ' ' + str(val)
        else:
            string = str(value)

        if self.Modify == True:
            for i in range(len(self.all_lines_raw)):
                line = self.all_lines_raw[i].strip()
                line = line.split(' ')
                if not check_block and 'BLOCK' in line and self.blockName in line:
                    check_block = True
                if not check_subBlock and 'SUBBLOCK' in line and self.subBlockName in line:
                    check_subBlock = True
                if key in line and check_block and check_subBlock:
                    check = True
                    index = self.all_lines_raw[i].index(key[0])
                    self.all_lines_raw[i] = self.all_lines_raw[i][
                                            :index] + key + ' = ' + string + '\n'  # WARNING: Change the " = " here if the separator between key and value is changed
                    break
        else:
            print('Function ReplaceFromKey can only be used in modify mode.')
            print('Please add "modify=True" while opening file.')
            print('Stopping.')
            sys.exit()
        if not check:
            print('WARNING : you tried to modify specifier ', key, ' but I did not find it.')
            if 'blockName' in self.__dict__.keys():
                print('You looked for Specifier ', key, ' in Block :', self.blockName)
            if 'subBlockName' in self.__dict__.keys():
                print('You looked for Specifier in SUBBLOCK :', self.subBlockName)
        return

    def Write(self):
        if self.Modify == True:
            self.file.seek(0)
            self.file.truncate()
            for line in self.all_lines_raw:
                self.file.write(line)
        else:
            print('Function Write can only be used in modify mode.')
            print('Please add "modify=True" while opening file.')
            print('Stopping.')
            sys.exit()
        return

    def ConvertLinesToFloat(self):
        check = True
        for i in range(len(self.currentLines)):
            try:
                self.currentLines[i] = self.currentLines[i].split(' ')
                if len(self.currentLines[i]) == 1:
                    self.currentLines[i] = self.currentLines[i][0]
                    self.currentLines[i] = float(self.currentLines[i])
                else:
                    for j in range(len(self.currentLines[i])):
                        self.currentLines[i][j] = float(self.currentLines[i][j])
            except:
                check = False
                self.currentLines[i] = ''
                pass
        if check == False:
            print('Conversion accomplished but it seems that for some lines it was not possible...')
        return

    def GetDictionaryFromBlock(self, KEYWORDS):
        dic = {}
        for key in self.GetKeysForBlock():
            typ = KEYWORDS.GetCharacterFromKey(key)
            if isinstance(typ, list):
                typ = typ[0]
            value = self.GetValueFromKey(key, typ)
            dic[key] = value
        return dic
