def callback(P):                                   #A function to validate the entry of numeric text boxes. Only allows number inputs. 
    if str.isdigit(P) or P == "":
        return True
    else:
        return False
    
def callbackFloat(P):                              #A function to validate the entry of numeric text boxes. Only allows float inputs.
    if P == "":
        return True
    try:
        float(P)
        return True
    except ValueError:
        return False