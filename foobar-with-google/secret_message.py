import base64

MESSAGE = '''
HEgXBg4MDhYARFJbQ0wAHQESGUhHRVQAHQ0PDgYIERZKT1FFVAYBFQYOCgoAVEFPTAAVBR0TFxhA T15TSgYFBgEGFggBBwJISFNKDggNGgYEBA4OCRtDU1dPTBAdDx0CCA4DSEhTSh0KBxEKBhJES11P QwAMCQ5CX0NVBwwEQE9eU0oYAgtSRA8=
'''

KEY = 'godsmokescrack'

result = []
for i, c in enumerate(base64.b64decode(MESSAGE)):
    #print(c)
    result.append(chr(c ^ ord(KEY[i % len(KEY)])))

print(''.join(result))
