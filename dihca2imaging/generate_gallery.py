#!/usr/bin/env python3
"""
Generate an HTML gallery for DIHCA spectral analysis diagnostic images.
This script automatically finds all diagnostic images in subdirectories and creates
a responsive web gallery for viewing them.
"""

import os
import glob
import json
from pathlib import Path

def find_diagnostic_images(results_directory):
    """Find all diagnostic PNG files in the results directory"""
    pattern = os.path.join(results_directory, "*", "*diagnostics.png")
    image_files = glob.glob(pattern)

    # Organize by source
    sources = {}
    for image_path in sorted(image_files):
        # Extract source name from directory
        rel_path = os.path.relpath(image_path, results_directory)
        source_dir = os.path.dirname(rel_path)

        if source_dir not in sources:
            sources[source_dir] = []

        sources[source_dir].append(rel_path.replace('\\', '/'))  # Ensure forward slashes for web

    return sources

def generate_gallery_html(results_directory, output_file=None):
    """Generate the HTML gallery file"""

    if output_file is None:
        output_file = os.path.join(results_directory, "gallery.html")

    # Find all diagnostic images
    sources = find_diagnostic_images(results_directory)

    if not sources:
        print("No diagnostic images found!")
        return

    # Generate JavaScript array for sources
    js_sources = []
    for source_name, images in sources.items():
        js_source = {
            "source": source_name,
            "images": images
        }
        js_sources.append(js_source)

    # Count totals
    total_sources = len(sources)
    total_images = sum(len(images) for images in sources.values())

    # HTML template
    html_content = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <meta content="IE=edge,chrome=1" http-equiv="X-UA-Compatible">
  <title>DIHCA Spectral Analysis Gallery</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 0;
      padding: 20px;
      background-color: #f5f5f5;
    }}
    .header {{
      margin: auto;
      max-width: 90%;
      text-align: center;
      background: white;
      padding: 20px;
      border-radius: 10px;
      margin-bottom: 20px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.1);
    }}
    .source-section {{
      margin: 20px 0;
      background: white;
      padding: 15px;
      border-radius: 10px;
      box-shadow: 0 2px 5px rgba(0,0,0,0.1);
    }}
    .source-title {{
      font-size: 1.5em;
      font-weight: bold;
      color: #333;
      margin-bottom: 15px;
      border-bottom: 2px solid #007acc;
      padding-bottom: 5px;
    }}
    .spw-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
      gap: 15px;
    }}
    .spw-item {{
      border: 1px solid #ddd;
      border-radius: 8px;
      overflow: hidden;
      background: #fafafa;
      transition: transform 0.2s, box-shadow 0.2s;
    }}
    .spw-item:hover {{
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
    }}
    .spw-title {{
      background: #007acc;
      color: white;
      padding: 8px 12px;
      font-weight: bold;
      font-size: 0.9em;
    }}
    .spw-image {{
      width: 100%;
      height: auto;
      cursor: pointer;
      transition: transform 0.2s;
      display: block;
    }}
    .spw-image:hover {{
      transform: scale(1.02);
    }}

    /* Modal styles */
    .modal {{
      display: none;
      position: fixed;
      z-index: 1000;
      left: 0;
      top: 0;
      width: 100%;
      height: 100%;
      background-color: rgba(0,0,0,0.9);
    }}
    .modal-content {{
      position: relative;
      margin: auto;
      padding: 0;
      width: 90%;
      max-width: 1200px;
      top: 50%;
      transform: translateY(-50%);
    }}
    .mySlides {{
      display: none;
      text-align: center;
    }}
    .mySlides img {{
      width: 100%;
      max-height: 80vh;
      object-fit: contain;
    }}
    .close {{
      position: absolute;
      top: 15px;
      right: 35px;
      color: #f1f1f1;
      font-size: 40px;
      font-weight: bold;
      cursor: pointer;
      z-index: 1001;
    }}
    .close:hover {{
      color: #bbb;
    }}
    .prev, .next {{
      cursor: pointer;
      position: absolute;
      top: 50%;
      width: auto;
      padding: 16px;
      margin-top: -50px;
      color: white;
      font-weight: bold;
      font-size: 20px;
      background: rgba(0,0,0,0.5);
      border: none;
      user-select: none;
      z-index: 1001;
    }}
    .next {{
      right: 0;
      border-radius: 3px 0 0 3px;
    }}
    .prev {{
      left: 0;
      border-radius: 0 3px 3px 0;
    }}
    .prev:hover, .next:hover {{
      background-color: rgba(0,0,0,0.8);
    }}
    .numbertext {{
      color: #f2f2f2;
      font-size: 12px;
      padding: 8px 12px;
      position: absolute;
      top: 0;
    }}
    .caption-container {{
      text-align: center;
      background-color: rgba(0,0,0,0.8);
      padding: 10px 16px;
      color: white;
      margin-top: 10px;
    }}
    .thumbnails {{
      display: flex;
      flex-wrap: wrap;
      justify-content: center;
      margin-top: 10px;
      max-height: 120px;
      overflow-y: auto;
    }}
    .thumbnail {{
      margin: 2px;
      cursor: pointer;
      opacity: 0.6;
      border: 2px solid transparent;
      border-radius: 4px;
    }}
    .thumbnail:hover, .thumbnail.active {{
      opacity: 1;
      border-color: #007acc;
    }}
    .thumbnail img {{
      width: 80px;
      height: 50px;
      object-fit: cover;
    }}

    .stats {{
      text-align: center;
      margin: 10px 0;
      color: #666;
      font-style: italic;
    }}

    .update-info {{
      text-align: center;
      margin: 10px 0;
      color: #888;
      font-size: 0.9em;
    }}
  </style>
</head>

<body>
  <div class="header">
    <h2>DIHCA Spectral Analysis Gallery</h2>
    <p>Diagnostic images for each source and spectral window. Click on any image to open the full-screen gallery browser, then use the left and right arrow keys to navigate between images.</p>
    <div class="stats">{total_sources} sources • {total_images} diagnostic images</div>
    <div class="update-info">Gallery generated automatically from spectral analysis results</div>
  </div>

  <div id="gallery-container">
    <!-- Content will be dynamically generated -->
  </div>

  <!-- Modal for full-screen viewing -->
  <div id="myModal" class="modal">
    <span class="close cursor" onclick="closeModal()">&times;</span>
    <div class="modal-content">
      <div id="slides-container"></div>
      <a class="prev" onclick="plusSlides(-1)">&#10094;</a>
      <a class="next" onclick="plusSlides(1)">&#10095;</a>
      <div class="caption-container">
        <p id="caption"></p>
      </div>
      <div class="thumbnails" id="thumbnails-container"></div>
    </div>
  </div>

  <script>
    let allImages = [];
    let slideIndex = 1;

    // Diagnostic images data (auto-generated)
    const diagnosticImages = {json.dumps(js_sources, indent=6)};

    // Generate the gallery content
    function generateGallery() {{
      const container = document.getElementById('gallery-container');
      let imageIndex = 0;

      diagnosticImages.forEach(sourceData => {{
        // Create source section
        const sourceSection = document.createElement('div');
        sourceSection.className = 'source-section';

        // Add source title
        const sourceTitle = document.createElement('div');
        sourceTitle.className = 'source-title';
        sourceTitle.textContent = sourceData.source;
        sourceSection.appendChild(sourceTitle);

        // Create grid for spectral windows
        const spwGrid = document.createElement('div');
        spwGrid.className = 'spw-grid';

        sourceData.images.forEach((imagePath, idx) => {{
          // Extract SPW info from filename
          const spwMatch = imagePath.match(/spw(\\d+)/);
          const spwNum = spwMatch ? spwMatch[1] : idx;

          // Create SPW item
          const spwItem = document.createElement('div');
          spwItem.className = 'spw-item';

          // Add SPW title
          const spwTitle = document.createElement('div');
          spwTitle.className = 'spw-title';
          spwTitle.textContent = `Spectral Window ${{spwNum}}`;
          spwItem.appendChild(spwTitle);

          // Add image
          const img = document.createElement('img');
          img.src = imagePath;
          img.className = 'spw-image';
          img.alt = `${{sourceData.source}} SPW ${{spwNum}}`;
          img.onclick = () => {{
            openModal();
            currentSlide(imageIndex + 1);
          }};

          // Handle image load errors
          img.onerror = function() {{
            this.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48cmVjdCB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgZmlsbD0iI2RkZCIvPjx0ZXh0IHg9IjUwJSIgeT0iNTAlIiBmb250LWZhbWlseT0iQXJpYWwiIGZvbnQtc2l6ZT0iMTQiIGZpbGw9IiM5OTkiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGR5PSIuM2VtIj5JbWFnZSBub3QgZm91bmQ8L3RleHQ+PC9zdmc+';
            this.alt = 'Image not found';
          }};

          spwItem.appendChild(img);
          spwGrid.appendChild(spwItem);

          // Add to global images array for modal
          allImages.push({{
            src: imagePath,
            alt: `${{sourceData.source}} SPW ${{spwNum}}`,
            caption: `${{sourceData.source}} - Spectral Window ${{spwNum}} Diagnostics`
          }});

          imageIndex++;
        }});

        sourceSection.appendChild(spwGrid);
        container.appendChild(sourceSection);
      }});

      // Generate modal content
      generateModalContent();
    }}

    // Generate modal content
    function generateModalContent() {{
      const slidesContainer = document.getElementById('slides-container');
      const thumbnailsContainer = document.getElementById('thumbnails-container');

      allImages.forEach((image, index) => {{
        // Create slide
        const slide = document.createElement('div');
        slide.className = 'mySlides';

        const numberText = document.createElement('div');
        numberText.className = 'numbertext';
        numberText.textContent = `${{index + 1}} / ${{allImages.length}}`;
        slide.appendChild(numberText);

        const img = document.createElement('img');
        img.src = image.src;
        img.alt = image.alt;
        slide.appendChild(img);

        slidesContainer.appendChild(slide);

        // Create thumbnail
        const thumbnail = document.createElement('div');
        thumbnail.className = 'thumbnail';
        thumbnail.onclick = () => currentSlide(index + 1);

        const thumbImg = document.createElement('img');
        thumbImg.src = image.src;
        thumbImg.alt = image.alt;
        thumbnail.appendChild(thumbImg);

        thumbnailsContainer.appendChild(thumbnail);
      }});
    }}

    // Modal functions
    function openModal() {{
      document.getElementById("myModal").style.display = "block";
      document.body.style.overflow = "hidden";
    }}

    function closeModal() {{
      document.getElementById("myModal").style.display = "none";
      document.body.style.overflow = "auto";
    }}

    function plusSlides(n) {{
      showSlides(slideIndex += n);
    }}

    function currentSlide(n) {{
      showSlides(slideIndex = n);
    }}

    function showSlides(n) {{
      const slides = document.getElementsByClassName("mySlides");
      const thumbnails = document.getElementsByClassName("thumbnail");
      const captionText = document.getElementById("caption");

      if (n > slides.length) {{ slideIndex = 1; }}
      if (n < 1) {{ slideIndex = slides.length; }}

      for (let i = 0; i < slides.length; i++) {{
        slides[i].style.display = "none";
      }}

      for (let i = 0; i < thumbnails.length; i++) {{
        thumbnails[i].classList.remove("active");
      }}

      if (slides.length > 0 && slideIndex > 0 && slideIndex <= slides.length) {{
        slides[slideIndex - 1].style.display = "block";
        thumbnails[slideIndex - 1].classList.add("active");
        captionText.innerHTML = allImages[slideIndex - 1].caption;
      }}
    }}

    // Keyboard navigation
    document.addEventListener('keydown', function(e) {{
      if (document.getElementById("myModal").style.display === "block") {{
        switch (e.key) {{
          case "ArrowLeft":
            plusSlides(-1);
            break;
          case "ArrowRight":
            plusSlides(1);
            break;
          case "Escape":
            closeModal();
            break;
        }}
      }}
    }});

    // Initialize gallery when page loads
    document.addEventListener('DOMContentLoaded', function() {{
      generateGallery();
      slideIndex = 1;
    }});
  </script>
</body>
</html>"""

    # Write the HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Gallery generated: {output_file}")
    print(f"Found {total_sources} sources with {total_images} diagnostic images")

    return output_file

def main():
    """Main function"""
    import sys

    if len(sys.argv) > 1:
        results_directory = sys.argv[1]
    else:
        results_directory = "spectral_analysis_results"

    if not os.path.exists(results_directory):
        print(f"Results directory not found: {results_directory}")
        return

    gallery_file = generate_gallery_html(results_directory)
    print(f"\\nOpen the gallery in your web browser:")
    print(f"file://{os.path.abspath(gallery_file)}")

if __name__ == "__main__":
    main()
