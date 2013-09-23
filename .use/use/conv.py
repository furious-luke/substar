import re, unicodedata
from xml.dom import minidom
from . import containers

try:
    import dateutil.parser
except:
    pass

# Don't assume googlemaps is available.
try:
    from googlemaps import GoogleMaps, GoogleMapsError
except:
    pass


def web_safe_translate(x):
    if x == ' ': return '-'
    if 'a' <= x <= 'z': return x
    if 'A' <= x <= 'Z': return x
    if '0' <= x <= '9': return x


def web_safe_make_translator():
    translations = ''.join(web_safe_translate(chr(c)) or chr(c) for c in range(256))
    deletions = ''.join(chr(c) for c in range(256) if web_safe_translate(chr(c)) is None)
    return translations, deletions


web_safe_translator = web_safe_make_translator()


def to_web_safe(value):
    if value is None:
        return None
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    return value.lower().translate(*web_safe_translator)


def strip_bad_unicode(value):
    if type(value) == type(u''):
        return value.replace(u'\u2019', '\'')
    return value


def to_iter(value):
    if hasattr(value, '__iter__') and not isinstance(value, dict):
        return value
    return [value]


def to_list(value):
    if value is None:
        return []
    elif isinstance(value, (str, bytes, dict)):
        return [value]
    elif hasattr(value, '__iter__'):
        return list(value)
    else:
        return [value]


def to_bool(value):
    if value is None:
        return None
    if isinstance(value, str):
        return not value.lower() in ['f', 'false']
    return bool(value)


def to_price(value, safe=False):
    if value is None:
        return None
    exp = r'\$\s*(\d+(?:\.\d+)?)'
    if isinstance(value, (str, bytes)):
        value = ' '.join(value.splitlines())
        if isinstance(value, bytes):
            match = re.search(exp, value, re.UNICODE)
        else:
            match = re.search(exp, value)
        if match:
            return float(match.group(1))
    try:
        return float(value)
    except:
        if safe:
            return None
        raise


def to_date(value):
    if value is None:
        return None
    dt = dateutil.parser.parse(value, fuzzy=True)
    return dt.date()


def to_time(value):
    if value is None:
        return None
    dt = dateutil.parser.parse(value, fuzzy=True)
    return dt.time()


def to_datetime(value):
    if value is None:
        return None
    dt = dateutil.parser.parse(str(value), fuzzy=True)
    return dt


def to_address(value, google_api_key, min_accuracy=1):
    if value is None:
        return None
    gmaps = GoogleMaps(google_api_key)
    if isinstance(value, (list, tuple)) and len(value) >= 2:
        if isinstance(value[0], (int, float)) and isinstance(value[1], (int, float)):
            result = gmaps.reverse_geocode(value[0], value[1])
        else:
            result = gmaps.geocode(value[0] + ' near ' + value[1])
    else:
        result = gmaps.geocode(value)
    address = result['Placemark'][0]['address']
    details = result['Placemark'][0]['AddressDetails']
    accuracy = details['Accuracy']
    # Raise an error if the accuracy is insufficient.
    if accuracy < min_accuracy:
        raise GoogleMapsError(602) # unknown address
    country = containers.find(details, 'CountryName', '')
    country_code = containers.find(details, 'CountryNameCode', '')
    state = containers.find(details, 'AdministrativeAreaName', '')
    state_code = containers.find(details, 'AdministrativeAreaNameCode', '')
    locality = containers.find(details, 'LocalityName', '')
    postal_code = containers.find(details, 'PostalCodeNumber', '')
    street_address = containers.find(details, 'ThoroughfareName', '')
    try:
        lng, lat = result['Placemark'][0]['Point']['coordinates'][0:2]
    except:
        lng, lat = None, None
    # Run through some common fixups.
    address_dict = dict(
        country=country,
        country_code=country_code,
        state=state,
        state_code=state_code,
        locality=locality,
        street_address=street_address,
        postal_code=postal_code,
        formatted=address,
        latitude=lat,
        longitude=lng,
    )
    address_fixups(address_dict)
    return address_dict


def address_fixups(address):
    if address['country'] == 'Australia':
        if address['state'] == 'VIC':
            address['state'] = 'Victoria'
            address['state_code'] = 'VIC'
        elif address['state'] == 'WA':
            address['state'] = 'Western Australia'
            address['state_code'] = 'WA'
        elif address['state'] == 'ACT':
            address['state'] = 'Australian Capital Territory'
            address['state_code'] = 'ACT'
        elif address['state'] == 'SA':
            address['state'] = 'South Australia'
            address['state_code'] = 'SA'
        elif address['state'] == 'QLD':
            address['state'] = 'Queensland'
            address['state_code'] = 'QLD'
        elif address['state'] == 'NSW':
            address['state'] = 'New South Wales'
            address['state_code'] = 'NSW'
        elif address['state'] == 'TAS':
            address['state'] = 'Tasmania'
            address['state_code'] = 'TAS'
        elif address['state'] == 'NT':
            address['state'] = 'Northern Territory'
            address['state_code'] = 'NT'


# def to_timezone(value, precise=False):
#     if value is None:
#         return None
#     latlng = to_latlng(value)
#     if precise:
#         url = "%s/%s/%s"%('http://www.earthtools.org/timezone', latlng[0], latlng[1])
#         response = urllib2.urlopen(url).read()
#         xml = minidom.parseString(response)
#         offset = float(xml.getElementsByTagName('offset')[0].childNodes[0].nodeValueue)
#     else:
#         offset = int(round(latlng[1]/15.0))
#     return offset


def to_url(value):
    # TODO: More checks/formatting.
    if value is None:
        return None
    return str(value)


def coerce_field(record, field, coer, **kwargs):
    if field in record:
        record[field] = coer(record[field], **kwargs)
    elif 'default' in kwargs:
        default = kwargs.pop('default')
        record[field] = coer(default, **kwargs)
    if field in record and record[field] is None:
        del record[field]


def int_field(record, field, **kwargs):
    coerce_field(record, field, int, **kwargs)


def unicode_field(record, field, **kwargs):
    coerce_field(record, field, unicode, **kwargs)


def bool_field(record, field, **kwargs):
    coerce_field(record, field, to_bool, **kwargs)


def price_field(record, field, **kwargs):
    coerce_field(record, field, to_price, **kwargs)


def date_field(record, field, **kwargs):
    coerce_field(record, field, to_date, **kwargs)


def time_field(record, field, **kwargs):
    coerce_field(record, field, to_time, **kwargs)


def datetime_field(record, field, **kwargs):
    coerce_field(record, field, to_datetime, **kwargs)


def address_field(record, field, **kwargs):
    coerce_field(record, field, to_address, **kwargs)


# def timezone_field(record, field, **kwargs):
#     coerce_field(record, field, to_timezone, **kwargs)


def url_field(record, field, **kwargs):
    coerce_field(record, field, to_url, **kwargs)
