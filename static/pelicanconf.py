AUTHOR = 'BdSlab'
SITENAME = 'BdSlab data'
SITEURL = ''

THEME = './themes/notmyidea'

PATH = 'content'

TIMEZONE = 'Europe/Brussels'

DEFAULT_LANG = 'en'
DEFAULT_ORPHANS = 10

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
#LINKS = (('Pelican', 'https://getpelican.com/'),
#         ('Python.org', 'https://www.python.org/'),
#         ('Jinja2', 'https://palletsprojects.com/p/jinja/'),
#         ('You can modify those links in your config file', '#'),)

LINKS = ()

# Social widget
#SOCIAL = (('You can add links in your config file', '#'),
#          ('Another social link', '#'),)
SOCIAL = ()
PAGINATED_TEMPLATES = {'index': None, 'tag': None, 'category': None, 'author': False}
DEFAULT_PAGINATION = False

# DIRECT_TEMPLATES = ['index']
# CATEGORY_SAVE_AS = ''

# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True