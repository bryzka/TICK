import configparser


def config(path="config"):
    c = configparser.ConfigParser()
    c.read(path)
    paths, opt = c["PATHS_req"], c["PATHS_opt"]
    return paths, opt

