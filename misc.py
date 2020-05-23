"""
Different miscellaneous procedures used in LONGRED package
"""

def message(text=None, noheader=False):
    """
    Display text message in adopted format (GUI field in a future)
    :param text: Message to display
    :param noheader: if True then name of routine will not be printer
    :return: None
    """
    header="LONGRED: "
    if noheader:
        header=''
    print("{}{}".format(header,text))