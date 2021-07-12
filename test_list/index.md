## Welcome to GitHub Pages

You can use the

## Welcome to GitHub Pages

You can use the
# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

{% assign pages = site.pages | sort: 'title' %}

<ul>
    {% for item in pages %}
        {% if item.dir == page.dir and item.path != page.path %}
           <li>  {{ item.name }} </li>
        {% endif %}
   {% endfor %}
</ul>

