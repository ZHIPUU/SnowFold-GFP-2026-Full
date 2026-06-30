from playwright.sync_api import sync_playwright

url='http://127.0.0.1:18082/'
with sync_playwright() as p:
    browser=p.chromium.launch(channel='msedge', headless=True)
    page=browser.new_page(viewport={"width":1440,"height":1000})
    page.goto(url, wait_until='networkidle')
    page.wait_for_timeout(3000)
    page.screenshot(path='D:/workspace/protein-design-atlas/artifacts/atlas_home.png', full_page=True)
    content = page.content()
    assert 'Protein Design Atlas' in content
    assert 'Sequences' in content
    page.get_by_text('序列库').click()
    page.wait_for_load_state('networkidle')
    page.wait_for_timeout(1000)
    assert 'Sequence Vault' in page.content()
    page.screenshot(path='D:/workspace/protein-design-atlas/artifacts/atlas_sequences.png', full_page=True)
    page.get_by_text('谱系网络').click()
    page.wait_for_timeout(2500)
    assert 'Topology Network' in page.content()
    page.screenshot(path='D:/workspace/protein-design-atlas/artifacts/atlas_network.png', full_page=True)
    browser.close()
print('PLAYWRIGHT_OK')
