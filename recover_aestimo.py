import json
import re

log_path = r'C:\Users\User\.gemini\antigravity\brain\b93d5945-3d95-4065-a321-5dc3a090b5da\.system_generated\logs\transcript_full.jsonl'
aestimo_path = r'C:\Users\User\Downloads\aestimo-feat-gui - simtest\aestimo.py'

with open(aestimo_path, 'r', encoding='utf-8') as f:
    content = f.read()

with open(log_path, 'r', encoding='utf-8') as f:
    for line in f:
        try:
            entry = json.loads(line)
        except:
            continue
        if 'tool_calls' in entry and entry['tool_calls']:
            for tc in entry['tool_calls']:
                args = tc.get('arguments', {})
                if args.get('TargetFile') == aestimo_path:
                    name = tc.get('name')
                    if name == 'replace_file_content':
                        target = args.get('TargetContent', '')
                        replacement = args.get('ReplacementContent', '')
                        if target in content:
                            content = content.replace(target, replacement)
                            print("Applied replace_file_content")
                    elif name == 'multi_replace_file_content':
                        chunks = args.get('ReplacementChunks', [])
                        for chunk in chunks:
                            target = chunk.get('TargetContent', '')
                            replacement = chunk.get('ReplacementContent', '')
                            if target in content:
                                content = content.replace(target, replacement)
                                print("Applied multi_replace_file_content chunk")

# Save recovered file
with open('aestimo_recovered.py', 'w', encoding='utf-8') as f:
    f.write(content)
print("Saved to aestimo_recovered.py")
