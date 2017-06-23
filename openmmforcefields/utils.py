from pkg_resources import resource_filename

def get_amber_directory():
    # TODO: Modify this to reflect where we install amber ffxml files
    return resource_filename('openmmforcefield', 'amber')

def get_charmm_directory():
    # TODO: Modify this to reflect where we install charmm ffxml files
    return resource_filename('openmmforcefield', 'charmm')
