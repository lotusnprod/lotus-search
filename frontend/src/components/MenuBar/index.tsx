'use client'

import Link from 'next/link';

import { useRouter } from 'next/router';
import { usePathname } from 'next/navigation'
const MenuBar: React.FC = () => {
  const pathname = usePathname()

  const menuItems = [
    { name: 'Home', path: '/' },
      { name: 'Structures', path: '/structures' },
      { name: 'Taxa', path: '/taxa' },
      { name: 'References', path: '/references' },
    // Add other menu items here
  ];

  return (
    <nav>
      <ul className="list-none flex justify-around">
        {menuItems.map((item) => (
          <li key={item.name} className={pathname === item.path ? "font-bold" : ""}>
            <Link href={item.path}>
              {item.name}
            </Link>
          </li>
        ))}
      </ul>
    </nav>
  );
};

export default MenuBar;